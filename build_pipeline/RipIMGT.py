"""
Rip IMGT/GENE-DB V/J germline sequences into local html + fasta files.

Source pages (GENElect query URLs), e.g.:
    https://www.imgt.org/genedb/GENElect?query=7.3+IGHV&species=Homo+sapiens

Output layout (unchanged -- FormatAlignments.py depends on it):
    IMGT_sequence_files/htmlfiles/<Species_underscore>_<GENETYPE>.html
    IMGT_sequence_files/fastafiles/<Species_underscore>_<GENETYPE>.fasta
e.g. Homo_sapiens_HV.fasta, Mus_BJ.fasta

IMGT now bot-blocks naive scraping (403/429) and the page HTML has changed
since 2015, so this script uses:
  - requests.Session with retries + exponential backoff on 403/429/5xx
  - User-Agent rotation per request
  - polite delay between requests (default 2 s)
  - tolerant FASTA extraction: BeautifulSoup <pre> blocks first,
    regex over the unescaped page text as fallback
  - resume/skip of already-downloaded non-empty files (use --force to redo)

Usage:
    python RipIMGT.py [--force] [--dry-run] [--delay 2.0]
                      [--only Homo_sapiens_HV,Mus_BJ] [--timeout 30]

Only stdlib + requests + bs4. No API keys.
"""

import argparse
import html as html_lib
import os
import random
import re
import sys
import time
from html.parser import HTMLParser
from html.entities import name2codepoint

try:
    import requests
    from requests.adapters import HTTPAdapter
    from urllib3.util.retry import Retry
except ImportError:  # pragma: no cover - caught at runtime with a clear message
    requests = None

try:
    from bs4 import BeautifulSoup
except ImportError:  # pragma: no cover
    BeautifulSoup = None


# Set globals
file_path = os.path.split(os.path.abspath(__file__))[0]
html_outpath = os.path.join(file_path, "IMGT_sequence_files", "htmlfiles")
fasta_outpath = os.path.join(file_path, "IMGT_sequence_files", "fastafiles")

# Define where to point the urls to.
# We have heavy, kappa, lambda, alpha, beta, gamma and delta chains.
# Both the v genes (imgt gapped amino acids) and the j genes (amino acids, are not gapped)

# Urls as of 04-12-14 (GENElect endpoint unchanged; page markup around it has)
urls = {"HV": "https://www.imgt.org/genedb/GENElect?query=7.3+IGHV&species=%s",
        "HJ": "https://www.imgt.org/genedb/GENElect?query=7.6+IGHJ&species=%s",
        "KV": "https://www.imgt.org/genedb/GENElect?query=7.3+IGKV&species=%s",
        "KJ": "https://www.imgt.org/genedb/GENElect?query=7.6+IGKJ&species=%s",
        "LV": "https://www.imgt.org/genedb/GENElect?query=7.3+IGLV&species=%s",
        "LJ": "https://www.imgt.org/genedb/GENElect?query=7.6+IGLJ&species=%s",
        "AV": "https://www.imgt.org/genedb/GENElect?query=7.3+TRAV&species=%s",
        "AJ": "https://www.imgt.org/genedb/GENElect?query=7.6+TRAJ&species=%s",
        "BV": "https://www.imgt.org/genedb/GENElect?query=7.3+TRBV&species=%s",
        "BJ": "https://www.imgt.org/genedb/GENElect?query=7.6+TRBJ&species=%s",
        "GV": "https://www.imgt.org/genedb/GENElect?query=7.3+TRGV&species=%s",
        "GJ": "https://www.imgt.org/genedb/GENElect?query=7.6+TRGJ&species=%s",
        "DV": "https://www.imgt.org/genedb/GENElect?query=7.3+TRDV&species=%s",
        "DJ": "https://www.imgt.org/genedb/GENElect?query=7.6+TRDJ&species=%s"

        # "HC": "https://www.imgt.org/genedb/GENElect?query=7.3+IGHC&species=%s",
        # "KC": "https://www.imgt.org/genedb/GENElect?query=7.3+IGKC&species=%s",
        # "LC": "https://www.imgt.org/genedb/GENElect?query=7.3+IGLC&species=%s",
        }


# Species as of 04-12-14
# Species as of 02-06-16 - alpaca added
# These are retrieved for all the antibodies
all_species = ["Homo+sapiens",
               "Mus",
               "Rattus+norvegicus",
               "Oryctolagus+cuniculus",
               "Macaca+mulatta",
               "Sus+scrofa",
               "Vicugna+pacos",
               "Bos+taurus"]


# These are retrieved for the tcr chains. There are a few more for gamma and delta chains but
# they are rare (structurally anyway) that it does not seem worth it.
all_tr_species = ["Homo+sapiens",
                  "Mus",
                  ]


# These do not have light chain sequences so we ignore (they're fish)
#           "Oncorhynchus+mykiss",
#           "Danio+rerio" ]

# Browser-like User-Agents, rotated per request so a single default UA
# doesn't get fingerprinted and 403-blocked immediately.
USER_AGENTS = [
    "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 "
    "(KHTML, like Gecko) Chrome/126.0.0.0 Safari/537.36",
    "Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:127.0) Gecko/20100101 Firefox/127.0",
    "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/605.1.15 "
    "(KHTML, like Gecko) Version/17.4 Safari/605.1.15",
    "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 "
    "(KHTML, like Gecko) Chrome/126.0.0.0 Safari/537.36",
]

# Matches one FASTA record in unescaped page text. Sequence lines may
# contain spaces/dots/gaps (IMGT gapped alignments); stop at next ">" or
# blank line. Header is everything after ">" up to end of line.
FASTA_RE = re.compile(
    r"^>(?P<header>[^\r\n]+)\r?$"
    r"(?P<body>(?:\r?\n(?![>])(?:[A-Za-z.\s*-]+))*)",
    re.MULTILINE,
)


def build_session(max_retries=5, backoff_factor=1.0):
    """requests.Session with retries + backoff on 403/429/5xx.

    403 is included because IMGT's bot protection raises it transiently;
    urllib3 honours Retry-After where the server sends one.
    """
    if requests is None:
        raise RuntimeError(
            "The 'requests' package is required. Install it with: pip install requests"
        )
    session = requests.Session()
    retry = Retry(
        total=max_retries,
        connect=max_retries,
        read=max_retries,
        status=max_retries,
        backoff_factor=backoff_factor,
        status_forcelist=[403, 429, 500, 502, 503, 504],
        allowed_methods=["GET"],
        respect_retry_after_header=True,
        raise_on_status=False,
    )
    adapter = HTTPAdapter(max_retries=retry)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    session.headers.update({
        "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
        "Accept-Language": "en-US,en;q=0.9",
        "Connection": "keep-alive",
    })
    return session


def parse_fasta_text(text):
    """Extract [(header, sequence)] records from plain (unescaped) text.

    Shared by the <pre>-block path and the whole-page regex fallback.
    Headers are returned WITHOUT the leading ">"; write_fasta adds it.
    """
    records = []
    name, chunks = None, []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None and chunks:
                records.append((name, "".join(chunks)))
            name = line[1:].strip()
            chunks = []
        elif name is not None:
            chunks.append(line.replace(" ", ""))
    if name is not None and chunks:
        records.append((name, "".join(chunks)))
    # Drop junk captures (IMGT headers are pipe-delimited accession lines)
    return [(h, s) for h, s in records if h and s and "|" in h]


def parse_sequences(htmlstring):
    """Tolerantly extract FASTA records from a GENElect HTML page.

    1. BeautifulSoup: scan <pre> blocks (where IMGT embeds the FASTA).
    2. Fallback: regex over the whole unescaped page text (covers markup
       changes where sequences are no longer inside <pre>).
    """
    unescaped = html_lib.unescape(htmlstring)

    if BeautifulSoup is not None:
        try:
            soup = BeautifulSoup(unescaped, "html.parser")
            records = []
            for pre in soup.find_all("pre"):
                records.extend(parse_fasta_text(pre.get_text()))
            if records:
                return records
        except Exception:
            pass  # fall through to regex fallback

    records = []
    for match in FASTA_RE.finditer(unescaped):
        header = match.group("header").strip()
        body = re.sub(r"\s+", "", match.group("body") or "")
        if header and body and "|" in header:
            records.append((header, body))
    if records:
        return records

    # Last resort: legacy line-based parse (handles pages with no "|" filter hit)
    records = []
    name, chunks = None, []
    for line in unescaped.splitlines():
        s = line.strip()
        if s.startswith(">"):
            if name is not None and chunks:
                records.append((name, "".join(chunks)))
            name, chunks = s[1:].strip(), []
        elif name is not None and s and re.fullmatch(r"[A-Za-z.\s*-]+", s):
            chunks.append(s.replace(" ", ""))
    if name is not None and chunks:
        records.append((name, "".join(chunks)))
    return records


# Html parser class (kept for backwards compatibility; parse_sequences
# above is the default path and handles post-2015 markup).
class GENEDBParser(HTMLParser):
    currenttag = None
    currentnamedent = None
    _data = []

    def handle_starttag(self, tag, attrs):
        self.currenttag = tag

    def handle_endtag(self, tag):
        self.currenttag = None

    def handle_data(self, data):
        split = data.split("\n")
        start = sum([1 if l[0] == ">" else 0 for l in split if len(l)])
        if self.currenttag == "pre" and (self.currentnamedent == ">" or start):
            # Two different ways of parsing the html based on how IMGT have formatted the pages.
            # For some reason they format gene db differently sometimes (legacy?)
            if start > 1:  # If you encounter more than one line in the data with a fasta ">" symbol, all sequences will be in the same packet
                name, sequence = None, ""
                for l in split:
                    if not l:
                        continue
                    if l[0] == ">":
                        if sequence:
                            self._data.append((name, sequence))
                            name, sequence = None, ""
                        name = l[1:].strip()
                    else:
                        sequence += l.replace(" ", "")
                if name and sequence:
                    self._data.append((name, sequence))
            else:  # Otherwise it will be done entry by entry
                try:
                    name = split[0]
                except IndexError:
                    return
                sequence = ("".join(split[1:])).replace(" ", "")
                if name.startswith(">"):
                    name = name[1:]
                self._data.append((name.strip(), sequence))

    def handle_entityref(self, name):
        self.currentnamedent = chr(name2codepoint[name])

    def handle_charref(self, name):
        if name.startswith('x'):
            self.currentnamedent = chr(int(name[1:], 16))
        else:
            self.currentnamedent = chr(int(name))

    def rip_sequences(self, htmlstring):
        """
        Method for this subclass that automates the return of data.
        Falls back to the tolerant parse_sequences on empty results.
        """
        self.reset()
        self._data = []
        self.currenttag = None
        self.currentnamedent = None
        self.feed(htmlstring)
        if not self._data:
            return parse_sequences(htmlstring)
        # Normalise legacy records (strip any leading ">" kept by old code)
        return [(n[1:].strip() if n.startswith(">") else n.strip(), s)
                for n, s in self._data if s]


parser = GENEDBParser()


def html_filename(species, gene_type):
    return os.path.join(html_outpath, "%s_%s.html" % (species.replace("+", "_"), gene_type))


def fasta_filename(species, gene_type):
    return os.path.join(fasta_outpath, "%s_%s.fasta" % (species.replace("+", "_"), gene_type))


def get_html(species, gene_type, force=False, session=None, timeout=30):
    """
    Get the html from IMGT. Returns the local filename, or False on failure.

    Skips download when a non-empty cached file exists unless force=True.
    403 -> likely bot-block (message says so); 429 -> rate-limited, back off.
    """
    filename = html_filename(species, gene_type)
    if os.path.isfile(filename) and os.path.getsize(filename) > 0 and not force:
        return filename
    if session is None:
        session = build_session()
    url = urls[gene_type] % species
    headers = {"User-Agent": random.choice(USER_AGENTS),
               "Referer": "https://www.imgt.org/"}
    try:
        response = session.get(url, headers=headers, timeout=timeout)
    except Exception as exc:
        print("Network error retrieving %s %s: %s" % (species, gene_type, exc),
              file=sys.stderr)
        return False
    if response.status_code == 403:
        print("Access denied (HTTP 403) for %s %s: IMGT bot protection blocked "
              "this client. Retry later, raise --delay, or download that page "
              "manually in a browser." % (species, gene_type), file=sys.stderr)
        return False
    if response.status_code == 429:
        print("Rate limited (HTTP 429) for %s %s: IMGT asked us to slow down. "
              "Increase --delay and retry." % (species, gene_type), file=sys.stderr)
        return False
    if response.status_code != 200 or not response.text:
        print("Bad response (HTTP %s) for %s %s" % (response.status_code, species, gene_type),
              file=sys.stderr)
        return False
    if ">" not in response.text or "GENElect" not in response.text and parse_sequences(response.text) == []:
        # No FASTA found: could be an empty query result (valid, e.g. rare
        # chain/species combo) or a changed/blocked page. Save anyway for
        # inspection but warn clearly.
        print("Warning: no FASTA records found in page for %s %s "
              "(empty result set or changed markup?)" % (species, gene_type),
              file=sys.stderr)
    os.makedirs(html_outpath, exist_ok=True)
    with open(filename, "w", encoding="utf-8", errors="replace") as outfile:
        outfile.write(response.text)
    return filename


def write_fasta(sequences, species, gene_type):
    """
    Write a fasta file containing all sequences.
    Same "<Species>_<GENETYPE>.fasta" layout FormatAlignments.py expects.
    """
    filename = fasta_filename(species, gene_type)
    os.makedirs(fasta_outpath, exist_ok=True)
    with open(filename, "w") as outfile:
        for name, sequence in sequences:
            print(">%s" % name, file=outfile)
            print(sequence, file=outfile)
    return filename


def ripfasta(species, gene_type, force=False, session=None, timeout=30):
    """
    Rip the fasta sequences for a species and gene type from IMGT.
    Returns 0 on success, 1 on failure. Skips existing non-empty fastas
    unless force=True (resume support).
    """
    fasta_file = fasta_filename(species, gene_type)
    if os.path.isfile(fasta_file) and os.path.getsize(fasta_file) > 0 and not force:
        return 0
    htmlfile = get_html(species, gene_type, force=force, session=session, timeout=timeout)
    if htmlfile:
        with open(htmlfile, encoding="utf-8", errors="replace") as infile:
            sequences = parse_sequences(infile.read())
        if not sequences:  # legacy parser as second opinion before giving up
            with open(htmlfile, encoding="utf-8", errors="replace") as infile:
                sequences = parser.rip_sequences(infile.read())
        if sequences:
            write_fasta(sequences, species, gene_type)
            return 0
        else:
            print("Bad parse: no sequences extracted for %s %s "
                  "(see %s)" % (species, gene_type, htmlfile), file=sys.stderr)
            return 1
    else:
        print("Bad Url: download failed for %s %s" % (species, gene_type), file=sys.stderr)
        return 1


def iter_targets(only=None):
    """Yield (species, gene_type) pairs, honouring skips and --only filter."""
    wanted = None
    if only:
        wanted = set()
        for item in only.split(","):
            item = item.strip().replace(" ", "")
            if not item:
                continue
            if "_" not in item:
                print("Ignoring malformed --only entry %r (want Species_GeneType, "
                      "e.g. Homo_sapiens_HV)" % item, file=sys.stderr)
                continue
            species, _, gene = item.rpartition("_")
            wanted.add((species.replace("_", "+"), gene))
    for gene_type in urls:
        for species in all_species:
            if gene_type[0] in "ABGD" and species not in all_tr_species:
                continue  # we don't want TCRs for all organisms
            if gene_type[0] in "KL" and species == "Vicugna+pacos":
                continue  # alpacas don't have light chains
            if wanted is not None and (species, gene_type) not in wanted:
                continue
            yield species, gene_type


def main(argv=None):
    """
    For all V and J gene types (H,K,L,A,B,G,D) parse IMGT database and extract fasta files.
    """
    argparser = argparse.ArgumentParser(description="Rip IMGT/GENE-DB germlines to fasta.")
    argparser.add_argument("--force", action="store_true",
                           help="Re-download and re-parse even if cached files exist.")
    argparser.add_argument("--dry-run", action="store_true",
                           help="List targets without downloading anything.")
    argparser.add_argument("--only", default=None,
                           help="Comma-separated subset, e.g. Homo_sapiens_HV,Mus_BJ.")
    argparser.add_argument("--delay", type=float, default=2.0,
                           help="Polite delay in seconds between requests (default 2.0).")
    argparser.add_argument("--timeout", type=float, default=30.0,
                           help="Per-request timeout in seconds (default 30).")
    args = argparser.parse_args(argv)

    targets = list(iter_targets(args.only))
    if args.dry_run:
        for species, gene_type in targets:
            print("%s %s -> %s" % (species, gene_type, fasta_filename(species, gene_type)))
        print("%d target(s), nothing downloaded (--dry-run)." % len(targets))
        return 0

    session = build_session()
    failures = 0
    for i, (species, gene_type) in enumerate(targets):
        if i > 0 and args.delay > 0:
            time.sleep(args.delay + random.uniform(0, 1.0))  # polite + jitter
        if ripfasta(species, gene_type, force=args.force,
                    session=session, timeout=args.timeout):
            print("Failed to retrieve %s %s" % (species, gene_type), file=sys.stderr)
            failures += 1
        else:
            print("Parsed and saved %s %s" % (species, gene_type))
    if failures:
        print("%d/%d target(s) failed." % (failures, len(targets)), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
