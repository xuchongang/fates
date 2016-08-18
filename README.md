# NGEE-T ED-CLM git repository
------------------------------

Repository for de-CLM-ification of ED. Contains a full copy of CLM
from recent svn trunk tags. Once ED has been isolated into a
standalone library, it will be moved into a new stand alone git repo
without a copy of CLM.

## Repository structure
-----------------------

Branches
  * master : ED-CLM trunk
  * clm-trunk : clm trunk tags from CESM svn
  * cesm2git : tools for pulling svn branch tags into git

## Returning changes to clm-svn
-------------------------------

Return changes to the CLM subversion repository by:

* determine the CLM trunk tag your version of the code is based on.
* Note that the script that pulls svn tags into git creates a separate
  commit and merge for squashed git externals, e.g. cime. This means
  that the commit id for a particular clm trunk tag may be one commit
  beyond the `pull clmX_Y_Z_rNNN` commit.
* determine the changeset id of the trunk version, and create a tag: clm-trunk-rXYZ
* determine the changeset id of the git version you want to return to subversion, and create a tag: ed-rABC
* create a patch file of the differences that should be moved to svn: `git diff clm-trunk-rXYZ ed-rABC > rXYZ-rABC.patch`
* create a new branch in subversion based on the SAME version of the trunk you used above.
* checkout a sandbox of your subversion branch.
* use the patch program to apply your patch to the subversion branch.
