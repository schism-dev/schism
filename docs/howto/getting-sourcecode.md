SCHISM modeling system is mostly distributed via github. To download the source code it is good idea to have `git` installed in your system. If you have git, you can download the sourcecode with the following commands - 

```bash
git clone https://github.com/schism-dev/schism.git
cd schism
```

This will give you the latest branch. If you want to get certain tagged versions use the `tag` command of git as `git checkout tags/v5.9.0` - for v5.9.0.

Note that if you are interested in the development, you might need to add your system’s ssh public key to your github.com account first before cloning. Alternatively, you can also go to https://github.com/schism-dev and directly download the zip files.

General users also have access to master branch of the development repository (https://github.com/schism-dev/schism), but this manual is currently for the tag version v5.9.0 only. An live documentation for the master version is currently in discussion. Currently, the web site and this manual is be updated when a newer version/bug fixes becomes available. 

Due to the large file size inside, the test suite used to test the code is still being distributed via svn. You’ll need svn v1.8 or above (see http://svnbook.red-bean.com/ for a manual for using svn). Svn clients on linux/unix/windows/Mac should all work. Following command will give you access to the verification tests - 

```bash
svn co https://columbia.vims.edu/schism/schism_verification_tests
```

Note that the test suite is used to test the master branch. Due to the differences between the release tags and master branch, you may need to modify some input files (e.g. param.nml) but we hope these tests  how you how we set up models for different problems.

Useful user info can be found in src/Readme.beta_notes (including change of format for input files and bug fixes).