from __future__ import print_function
from builtins import next
import shutil
from distutils.core import setup, Command
import os
import glob
import errno
import subprocess


__author__ = 'pyrocko devs'

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


def check_broken_links(filenames, undangle):
    for fn in filenames:
        try:
            os.stat(fn)
        except OSError as e:
            if e.errno == errno.ENOENT:
                if not undangle:
                    print(
                        'file %s does not exist or is a broken symlink'
                        % fn)
                    print(
                        'broken symlinks can be removed using --undangle')
                else:
                    os.unlink(fn)
                    print('Unlinked file:  %s' % fn)


class PassSetup(SetupBuildCommand):
    descrition = """install doesnt work. Use "python setup.py link" instead."""
    user_options = []

    def run(self):
        print("""install doesnt work. Use "python setup.py link" instead.""")


class LinkSnufflingFiles(SetupBuildCommand):

    description = "Create symbolic links and subdirectory in $HOME/.snufflings"

    user_options = [
        ('force', None, 'force overwriting of existing symbolic links.'),
        ('choice=', None, 'Comma separated list of snufflings to link')
    ]

    def initialize_options(self):
        self.force = False
        self.undangle = True
        self.choice = []
        self.excluded_dirs = ['.git', 'screenshots']

    def get_target_files(self):
        cwd = os.getcwd()
        files = []
        if self.choice:
            choices = self.choice.split(',')
            files = []
            for c in choices:
                files.extend(glob.glob(pjoin(cwd, c)))

        else:
            files = glob.glob(pjoin(cwd, '*.py'))
            subs = next(os.walk('.'))[1]
            subs = [x for x in subs if x[0] != '.']
            for excl in self.excluded_dirs:
                if excl in subs:
                    subs.remove(excl)
            for sub in subs:
                files.append(pjoin(cwd, sub))
            files.remove(cwd+'/setup.py')
        return files

    def build_extensions(self, root_dir):
        for roots, dirs, _files in os.walk(root_dir, topdown=True):
            dirs[:] = [d for d in dirs if d not in self.excluded_dirs]

            # okada temporarily disabled: causes segfaults on ubuntu 14.04 when
            # snuffler is closed.
            if 'Makefile' in _files and not 'okada.py' in _files:
                print('\nbuilding %s' % roots)
                try:
                    subprocess.check_call(['make', '-C', roots])
                except subprocess.CalledProcessError:
                    print(' failed building %s' % (roots))

    def run(self):
        cwd = os.getcwd()
        snufflings = pjoin(os.getenv('HOME'), '.snufflings')

        cache = pjoin(snufflings, '__pycache__')
        if os.path.exists(cache):
            shutil.rmtree(cache)

        # look for dangling symlinks inside .snufflings:
        check_broken_links(glob.glob(snufflings+'/*'), self.undangle)
        self.build_extensions(cwd)
        files = self.get_target_files()
        for fn in files:
            try:
                target = fn.replace(cwd, snufflings)
                os.symlink(fn, target)
            except OSError as e:
                if e.errno == errno.EEXIST:
                    if os.path.islink(target):
                        if self.force:
                            os.remove(target)
                            os.symlink(fn, target)
                        else:
                            print('target exists: ', target)
                    else:
                        print('target exists and is not a symbolic link: ',
                              target)


setup(name='contrib-snufflings',
      version='1.0',
      description='User contributed snufflings.',
      packages=[],
      cmdclass={'link': LinkSnufflingFiles,
                'install': PassSetup}
      )
