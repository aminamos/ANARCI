import shutil, os, subprocess
import site, sys
from setuptools import setup
from setuptools.command.install import install

PKG_DIR = os.path.join('lib', 'python', 'anarci')

class CustomInstallCommand(install):
   def run(self):
       install.run(self)
       # Post-installation routine
       ANARCI_LOC = os.path.join(site.getsitepackages()[0], 'anarci') # site-packages/ folder
       ANARCI_BIN = sys.executable.split('python')[0] # bin/ folder

       shutil.copy('bin/ANARCI', ANARCI_BIN) # copy ANARCI executable
       print("INFO: ANARCI lives in: ", ANARCI_LOC) 

       # Prebuilt germlines.py and HMMs ship with the package (built with
       # pyhmmer, no HMMER binary required). Only rebuild from IMGT when they
       # are missing, e.g. in a bare source checkout without the data.
       prebuilt = (os.path.exists(os.path.join(PKG_DIR, 'germlines.py')) and
                   os.path.exists(os.path.join(PKG_DIR, 'dat', 'HMMs', 'ALL.hmm')))
       if prebuilt:
           print('INFO: using prebuilt germlines and HMMs.')
       else:
           # Build HMMs from IMGT germlines
           os.chdir("build_pipeline")
           print('INFO: Downloading germlines from IMGT and building HMMs...')
           print("INFO: running 'RUN_pipeline.sh', this will take a couple a minutes.")
           proc = subprocess.Popen(["bash", "RUN_pipeline.sh"], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
           o, e = proc.communicate()

           print(o.decode())
           print(e.decode())
           
           # Copy HMMs where ANARCI can find them
           shutil.copy( "curated_alignments/germlines.py", ANARCI_LOC )
           dat = os.path.join(ANARCI_LOC, "dat")
           if not os.path.isdir(dat):
               os.mkdir(dat)
           shutil.copytree( "HMMs", os.path.join(ANARCI_LOC, "dat", "HMMs"), dirs_exist_ok=True )
          
           # Remove data from HMMs generation
           try:
               shutil.rmtree("curated_alignments/")
               shutil.rmtree("muscle_alignments/")
               shutil.rmtree("HMMs/")
               shutil.rmtree("IMGT_sequence_files/")
           except OSError:
               pass

setup(name='anarci',
     version='1.3',
     description='Antibody Numbering and Receptor ClassIfication',
     author='James Dunbar',
     author_email='opig@stats.ox.ac.uk',
     url='http://opig.stats.ox.ac.uk/webapps/ANARCI',
     packages=['anarci'],
     package_dir={'anarci': 'lib/python/anarci'},
     package_data={'anarci': ['dat/HMMs/*']},
     data_files = [ ('bin', ['bin/muscle', 'bin/muscle_macOS', 'bin/ANARCI']) ],
     include_package_data = True,
     install_requires=['pyhmmer'],
     scripts=['bin/ANARCI'],
     cmdclass={"install": CustomInstallCommand, }, # Run post-installation routine
    )
