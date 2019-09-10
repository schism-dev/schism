# Contributing developers: please read

## SCHISM github Workflow (Trunk Based Development)

The following is tested on linux git. Should work same on other platforms.

### Clone
In git, there are multiple repositories (repos), one local and multiple ‘remotes’. Different git commands are used for these. The local is your private repo, and until you publish your changes there (e.g. via ‘push’), others cannot access it.

Also note that most git commands can be issued anywhere under the top dir as they operate on the entire repo.

To clone a copy of the repository as an external user - simply do the following

```bash
$ mkdir schism_git #or whatever name you want
$ cd schism_git
$ git clone clone https://github.com/schism-dev/schism.git
$ cd schism
```

This will not ask for any username-password combo. If you would like to push, you have to provide your username and password to validate your write access to the repository.

For password-less ssh-based authentication, you need to first add ssh key (e.g. .ssh/id_rsa.pub) to your github account under ‘Settings->SSH and GPG keys’.

```bash
$ mkdir schism_git #or whatever name you want
$ cd schism_git
$  git clone https://github.com/schism-dev/schism.git
# this clones the official repo, and points local branch to origin/master branch which is like svn’s trunk
$  cd schism
$ git remote set-url origin git@github.com:schism-dev/schism.git
```
We suggest to use ssh instead of https form to get around issues with https. Which protocol is supported depends on your machine and firewall configuration.
```
$ git remote -v  #examine remotes

origin  git@github.com:schism-dev/schism.git (fetch)
origin  git@github.com:schism-dev/schism.git (push)
```
Examine all branches, local and remote using `git branch -a`. Note the pointer ‘HEAD’ now points to origin/master. There are other branches in origin (remote) that maybe created by other developers)
```bash
$ git branch -a

  * master
  remotes/origin/HEAD -> origin/master
  remotes/origin/cmake
  remotes/origin/master
  remotes/origin/new_branch
```
### Working with repo
After cloning as explained above your working directory (local HEAD) will be pointing to origin/master.
In principle, you can create changes and push directly to the master. 

Option 1: Direct push
This option is good for small incremental/non-controversial  changes.
(edit local files and then commit into local repo first (not remote!))

```bash
$ cd src/Hydro; vi schism_newmodule.F90 #make some changes
$ git add schism_newmodule.F90 #since this is a new file, need to add first
$ git commit –av  
# ‘-a’ bypasses staging step
# this cmd can be issued in any dir under schism/ and it will actually commit all new changes to local repo
$ git pull #pull others’ changes from official repo; do this often to be in sync with remote
$ git status # useful for examining uncommitted files etc. 
```
Note that .gitignore is a great way to ignore some intermediate files such as *.o etc)
When ready to publish your changes from local repo to remote, make sure your working copy is clean, i.e. nothing to commit into local repo. Then push)

```bash
$ git push origin master
```

Option 2: branching inside official repo
Suppose you want to collaborate with another colleague on some complex code changes that would break the master if you pushed them directly. For that reason we commit intermediary changes that break to code to branches (also called “feature branches”). In the case above, suppose you want to create a new  branch.

```
$ pwd

schism
```
```
$ git checkout master #starting from master branch
$ git pull #update from repo
$ git branch new_branch master #create a new branch ‘new_branch’ based on master
$ git checkout new_branch #switch to new_branch
$ git status #examine)

# On branch new_branch
nothing to commit (working directory clean)
(Push to remote to establish your new branch in remote repo so others can access it)

`$ git push --set-upstream origin new_branch`

(now do similar things as in Option 1, edit and commit to this branch locally. Push to remote when ready)

`$ git push` (publish your changes to the remote origin/new_branch)

When ready to merge your branch with master, there are two approaches: via pull request (PR) or local merge+push

3.1 Local merge + push

```bash
$ git checkout master #switch to master branch
$ git pull origin master   #or:  git reset --hard origin/master  (bring it up to date)
$ git merge --squash new_branch #merge ‘new_branch’ to master; resolve conflict if necessary; write a single commit message
$ git push origin master #push to master branch of remote repo

Clean up: 
(need to switch to another branch than ‘new_branch’’ first by using git checkout...)

```bash
$ git branch -D new_branch #delete from local repo
$ git push origin --delete new_branch #remove from remote repo
```
To make the commit log clean, a useful way is to rebase (i.e. move HEAD around) so that commits into new_branch are squashed in the log. The best way to rebase is to use github’s PR.

3.2 Pull request
On github, submit a PR, so other developers can examine the changes and provide feedback. When consensus is reached, some developer can do final merge on github as well, with optional rebase. The GUI there provides good help on all of these tasks.

Option 3: fork
