from distutils.file_util import copy_file
from distutils.core import setup, Command
import os
import glob
import errno

__author__ = 'marius'

pjoin = os.path.join


class SetupBuildCommand(Command):
    """
    Master setup build command to subclass from.
    """
    def initialize_options(self):
        """
        Setup the current dir.
        """
        self._dir = os.getcwd()

    def finalize_options(self):
        """
        Set final values for all the options that this command supports.
        """


class PassSetup(SetupBuildCommand):
    descrition = """install doesnt work. Use "python setup.py link" instead.""" 
    user_options = []
    def run(self):
        print """install doesnt work. Use "python setup.py link" instead."""

class LinkSnufflingFiles(SetupBuildCommand):

    description = "Create symbolic links and subdirectory in $HOME/.snufflings"

    user_options = [('force', None, 
                            'force overwriting of existing symbolic links.'),
                    ('choice=', None, 'Comma separated list of snufflings to link'), ]

    def initialize_options(self):
        self.force = False
        self.choice = []

    def run(self):
        snufflings = pjoin(os.getenv('HOME'), '.snufflings')
        cwd = os.getcwd()
        if self.choice:
            choices = self.choice.split(',')
            files = []
            for c in choices:
                files.extend(glob.glob(pjoin(cwd, c)))

        else:
            files = glob.glob(pjoin(cwd, '*'))

        for fn in files:
            try:
                target = fn.replace(cwd, snufflings)
                os.symlink(fn, target)
            except OSError, e:
                if e.errno==errno.EEXIST:
                    if os.path.islink(target):
                        if self.force:
                            os.remove(target)
                            os.symlink(fn, target)
                        else:
                            print 'target exists: ', target
                    else:
                        print 'target exists and is not a symbolic link: ',target


setup(name='contrib-snufflings',
      version='1.0',
      description='User contributed snufflings.',
      packages=[],
      cmdclass={'link': LinkSnufflingFiles,
                'install': PassSetup}
      )
