SCHISM modeling system is mostly distributed via github. To download the source code it is good idea to have `git` installed in your system. If you have git, you can download the sourcecode with the following commands - 

```bash
git clone https://github.com/schism-dev/schism.git
cd schism
```

This will give you the latest master branch. If you want to get certain tagged versions use the `tag` command of git as (for example) `git checkout tags/v5.10.0` - for v5.10.0.

Note that sometimes you might need to add your system’s ssh public key to your github.com account first before cloning. 
If cloning fails, you can also go to https://github.com/schism-dev and directly download the zip files.

General users have access to all branches and tags  of the repository (https://github.com/schism-dev/schism). You can find online manuals for the latest stable
 master or newer tags after v5.9.0. The developers will constantly update the manuals to keep them up to date as much as possible. 

Due to the large file size inside, the test suite used to test the code is still being distributed via svn. You’ll need svn v1.8 or above (see http://svnbook.red-bean.com/ for a manual for using svn). Svn clients on linux/unix/windows/Mac should all work. Following command will give you access to the verification tests - 

```bash
svn co https://columbia.vims.edu/schism/schism_verification_tests
```

Note that the test suite is used to test the latest master branch. Due to the differences 
between the release tags and master branch, you may need to modify some input files (e.g. `param.nml`) 
but we hope these tests  how you how we set up models for different problems.

Useful user info can be found in src/Readme.beta_notes (including change of format for input files and bug fixes).
