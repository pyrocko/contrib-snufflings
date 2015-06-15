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
                    ('choice=', None, 'Comma separated list of snufflings to link'),
                    ('undangle', None, 'Unlink broken (dangling) symlinks in $HOME/.snufflings'), ]

    def initialize_options(self):
        self.force = False
        self.undangle = False
        self.choice = []

    def run(self):
        snufflings = pjoin(os.getenv('HOME'), '.snufflings')
        cwd = os.getcwd()

        # look for dangling symlinks inside .snufflings:
        for fn in glob.glob(snufflings+'/*'):
            try:
                os.stat(fn)
            except OSError, e:
                if e.errno == errno.ENOENT:
                    if not self.undangle:
                        print 'file %s does not exist or is a broken symlink' % fn
                        print 'broken symlinks can be removed using --undangle'
                    else:
                        os.unlink(fn)
                        print 'Unlinked file:  %s'%fn

        if self.choice:
            choices = self.choice.split(',')
            files = []
            for c in choices:
                files.extend(glob.glob(pjoin(cwd, c)))

        else:
            files = glob.glob(pjoin(cwd, '*.py'))
            subs = next(os.walk('.'))[1]
            subs = filter(lambda x: x[0]!='.', subs)
            subs.remove('screenshots')
            for sub in subs:
                files.append(pjoin(cwd, sub))
            files.remove(cwd+'/setup.py')

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
