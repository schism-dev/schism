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

This will not ask for any username-password combo. If you would like to push, you need to have developers credentials and have to provide your username and password to validate your write access to the repository.

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
Examine all branches, local and remote using `git branch -a`. Note the pointer ‘HEAD’ now points to origin/master. There are other branches in origin (remote) that may be created by other developers)
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

Since git uses distributed repos, ‘pull’ from remote (i.e. fetch and merge with your local copy) will incorporate changes from others into your local repo, which is considered a ‘delta’ that will be recorded in commit history log when you try to ‘push’ back into remote. This is burdensome for reviewers because the merge does not represent your own commits. One way to squash this ‘delta’ is to rebase, which carries its own risk of entangling commit history and thus divergence when multiple users are working on a same branch. Therefore, as a general rule, if git asks you for a commit message during pull, please indicate it as ‘PULL MERGE’ so reviewers can safely ignore.

Option 1: Direct push

This option is good for small incremental/non-controversial  changes.

(edit local files and then commit into local repo first (not remote!))

```bash
$ cd src/Hydro; vi schism_newmodule.F90 #make some changes
$ git add schism_newmodule.F90 #since this is a new file, need to add first
$ git commit –av  # ‘-a’ bypasses staging step
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
$ git checkout master #starting from master branch
$ git pull #update from repo
$ git branch new_branch master #create a new branch ‘new_branch’ based on master
$ git checkout new_branch #switch to new_branch
$ git status #examine
# On branch new_branch
nothing to commit (working directory clean)
```

(Push to remote to establish your new branch in remote repo so others can access it)

`$ git push --set-upstream origin new_branch`

(now do similar things as in Option 1, edit and commit to this branch locally. Push to remote when ready)

`$ git push --set-upstream origin new_branch`

(Get updates from other with pull)
`$ git pull origin new_branch`

When ready to merge your branch with master, there are two approaches: via pull request (PR) or local merge+push

2.1 Local merge + push

```bash
$ git checkout master #switch to master branch
$ git pull origin master   #or:  git reset --hard origin/master  (bring it up to date)
$ git merge --squash new_branch #merge ‘new_branch’ to master; resolve conflict if necessary; write a single commit message
$ git commit -av 
$ git push origin master #push to master branch of remote repo
```

Clean up: 
(need to switch to another branch than ‘new_branch’’ first by using git checkout...)

```bash
$ git branch -D new_branch #delete from local repo
$ git push origin --delete new_branch #remove from remote repo
```
To make the commit log clean, a useful way is to rebase (i.e. move HEAD around) so that commits into new_branch are squashed in the log. The best way to rebase is to use github’s PR.

2.2 Pull request

On github, submit a PR, so other developers can examine the changes and provide feedback. When consensus is reached, some developer can do final merge on github as well, with optional rebase. The GUI there provides good help on all of these tasks.

Option 3: fork

If you intend to keep your branch for a while, you may consider forking off from the official repo.
However, beware that this can lead to divergence over time unless you diligently try to merge.
To fork, go to github.com and use the fork function there. The fork will have your name in the new repo.
Use PR function to manage merge.

### Other useful commands
#### To checkout a specific commit
```bash
1) Need to know the commit number (SHA1 aka hash); first 7 digits are sufficient (e.g. 40b5ad0cbd26d026caf934bff9c12723e7773f65)
2) git clone into a separate dir. Better rename the dir to some meaningful name
3) git describe 40b5ad0 #(outputs a SHA1 that can be used to checkout)
r5255-43-g40b5ad0
4) git checkout r5255-43-g40b5ad0 
```

#### Misc
```bash
$ git log --pretty=format:"%h %s" --graph --all #(useful for examining branch)
$ git log -p #(show all changes in all files)
$ git log --since=2.weeks
$ git rm -f ... and then commit #(remove a file/dir; push will remove working copy as well as repository copy)
$ git mv f1 f2 #rename
$ git checkout -- <files> #(like svn revert (overwrite working copy with repo copy; <files> can be a dir)
```
