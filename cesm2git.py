#!/usr/bin/env python
"""create "shallow" git clones of cesm by pulling in specified svn
branch tags.

Author: Ben Andre <andre@ucar.edu>

"""

from __future__ import print_function

import sys

if sys.hexversion < 0x02070000:
    print(70 * "*")
    print("ERROR: {0} requires python >= 2.7.x. ".format(sys.argv[0]))
    print("It appears that you are running python {0}".format(
        ".".join(str(x) for x in sys.version_info[0:3])))
    print(70 * "*")
    sys.exit(1)

import argparse
import os
import shutil
import subprocess
import time
import traceback

if sys.version_info[0] == 2:
    from ConfigParser import SafeConfigParser as config_parser
else:
    from configparser import ConfigParser as config_parser


# -------------------------------------------------------------------------------
#
# User input
#
# -------------------------------------------------------------------------------

def commandline_options():
    """Process the command line arguments.

    """
    parser = argparse.ArgumentParser(
        description='FIXME: python program template.')

    parser.add_argument('--backtrace', action='store_true',
                        help='show exception backtraces as extra debugging '
                        'output')

    parser.add_argument('--debug', action='store_true',
                        help='extra debugging output')

    parser.add_argument('--config', nargs=1, required=True,
                        help='path to config file')

    parser.add_argument('--repo', nargs=1, default=['ed-clm-git'],
                        help='path to ed-clm-git repo, relative to cwd.')

    parser.add_argument('--feelin-lucky', action='store_true', default=False,
                        help='push the updated changes back to the master repo')

    options = parser.parse_args()
    return options


def read_config_file(filename):
    """The configuration file contains information about which cesm tag to
    check out and which svn externals need to be changed:

    [cesm]
    repo = https://svn-ccsm-models.cgd.ucar.edu
    piorepo = https://parallelio.googlecode.com
    tag = cesm1/alphas/tags/cesm1_3_alpha13b

    [externals]
    scripts = scripts/trunk_tags/scripts4_140915
    models/ocn/pop2  = pop2/trunk_tags/cesm_pop_2_1_20140828



    NOTE: the options in the external sections are exactly the
    directory paths used in SVN_EXTERNAL_DIRECTORIES, the values are
    the url minus the repo prefix.

    FIXME(bja, 201410) This assumes we won't want to update an
    external that isn't stored in svn. Probably a bad assumption.

    """
    print("Reading configuration file : {0}".format(filename))

    cfg_file = os.path.abspath(filename)
    if not os.path.isfile(cfg_file):
        raise RuntimeError("Could not find config file: {0}".format(cfg_file))

    config = config_parser()
    config.read(cfg_file)

    repo_config = {}

    def _check_for_required_section(conf, section):
        if not conf.has_section(section):
            raise RuntimeError(
                "ERROR: repo config file must contain a '{0}' section".format(section))

    def _get_section_required_option(conf, section, option, upper_case=False):
        sect = list_to_dict(conf.items(section))
        if option not in sect and option.upper() not in sect:
            raise RuntimeError("ERROR: repo config section '{0}' must contain"
                               "a '{1}' keyword.".format(section, option))

        opt = {option: sect[option]}
        if upper_case is True:
            opt = {option.upper(): sect[option]}

        return opt

    section = "git"
    repo_config[section] = {}
    _check_for_required_section(config, section)
    keys = ["branch", ]
    for k in keys:
        repo_config[section].update(
            _get_section_required_option(
                config, section, k))

    section = "cesm"
    repo_config[section] = {}
    _check_for_required_section(config, section)
    keys = ["repo", "tag"]
    for k in keys:
        repo_config[section].update(
            _get_section_required_option(
                config,
                section,
                k))

    section = "externals"
    repo_config[section] = {}
    if config.has_section(section):
        repo_config[section].update(list_to_dict(config.items(section)))

    return repo_config


# -------------------------------------------------------------------------------
#
# misc work functions
#
# -------------------------------------------------------------------------------
def list_to_dict(input_list, upper_case=False):
    output_dict = {}
    for item in input_list:
        key = item[0]
        value = item[1]
        if upper_case is True:
            key = key.upper()
        output_dict[key] = value
    return output_dict


def new_tag_from_config(config):
    """Generate a meaningful, if very verbose tag name from the specified
    cesm tag and externals

    """
    new_tag = config["cesm"]["tag"].split('/')[-1]
    for ext in config["externals"]:
        tag = config["externals"][ext].split('/')[-1]
        new_tag += "-{0}".format(tag)

    print("Creating new tag: {0}".format(new_tag))
    return new_tag


def remove_current_working_copy():
    """Removes the current working copy of cesm so that svn checkout will work.
    """

    rm_files = [
        "ChangeLog",
        "ChangeLog_template",
        "Copyright",
        "README",
        "SVN_EXTERNAL_DIRECTORIES",
        "ExpectedTestFails.xml",
    ]
    for f in rm_files:
        if os.path.exists(f):
            os.remove(f)

    rm_dirs = [
        "components",
        "cime",
    ]

    for d in rm_dirs:
        if os.path.exists(d):
            shutil.rmtree(d)


# -------------------------------------------------------------------------------
#
# svn wrapper functions
#
# -------------------------------------------------------------------------------
def svn_checkout_cesm(cesm_config, debug):
    """Checkout the user specified cesm tag
    """
    print("Checking out cesm tag from svn...", end='')
    url = cesm_config['repo']
    cesm_tag = cesm_config['tag']
    cmd = [
        "svn",
        "co",
        "{0}/{1}".format(url, cesm_tag),
        "."
    ]
    output = subprocess.STDOUT
    if debug:
        print("\n")
        print(" ".join(cmd))
        output = None
    subprocess.check_output(cmd, shell=False, stderr=output)
    if not debug:
        print(" done.")


def update_svn_externals(temp_repo_dir, repo_url, external_mods):
    """Backup the svn externals file, read it in and modify according to
    the user config, then write the new externals file.

    """
    print("Updating svn externals...", end='')
    externals_filename = "{0}/SVN_EXTERNAL_DIRECTORIES".format(temp_repo_dir)
    shutil.copy2(externals_filename, "{0}.orig".format(externals_filename))

    new_externals = []
    with open(externals_filename, 'r') as externals_file:
        externals = externals_file.readlines()
    for line in externals:
        temp = line.split()
        if len(temp) == 2:
            ext = temp[0]
            for e in external_mods:
                if e.strip() == ext:
                    line = "{0}            {1}/{2}\n".format(
                        ext,
                        repo_url,
                        external_mods[e])
        new_externals.append(line)

    with open(externals_filename, 'w') as externals_file:
        for line in new_externals:
            externals_file.write(line)

    svn_set_new_externals()
    svn_update("components")
    print(" done.")


def svn_set_new_externals():
    """
    """
    cmd = [
        "svn",
        "propset",
        "svn:externals",
        "-F",
        "SVN_EXTERNAL_DIRECTORIES",
        ".",
    ]
    subprocess.check_output(cmd, shell=False, stderr=subprocess.STDOUT)


def svn_update(path):
    """
    """
    cmd = [
        "svn",
        "update",
        path,
    ]
    subprocess.check_output(cmd, shell=False, stderr=subprocess.STDOUT)


def svn_switch(temp_repo_dir, switch_dir, url, tag):
    """run svn switch to update an external to a different tag

    FIXME(bja, 20141008): this means the SVN_EXTERNAL_DIRECTORIES are
    wrong. Need to update the externals file then update so everything
    is in sync!

    """
    os.chdir("{0}/{1}".format(temp_repo_dir, switch_dir))
    cmd = [
        "svn",
        "switch",
        "{0}/{1}".format(url, tag),
    ]
    subprocess.check_output(cmd, shell=False, stderr=subprocess.STDOUT)
    os.chdir(temp_repo_dir)


# -------------------------------------------------------------------------------
#
# git wrapper functions
#
# -------------------------------------------------------------------------------
def clone_cesm_git(repo_dir, temp_repo_dir):
    """Clone the existing git repo.

    NOTE: assumes a fixed directory structure. If this script is
    executed from directory 'work', then work/ed-clm-git is the main git
    repo to pull new src into.

    """
    print("Cloning git repo at : {0}".format(repo_dir))
    cmd = [
        "git",
        "clone",
        repo_dir,
        temp_repo_dir,
    ]
    subprocess.check_output(cmd, shell=False, stderr=subprocess.STDOUT)


def switch_git_branch(branch):
    """All cesm changes from upstream svn are pulled onto the git 'trunk' branch
    """
    cmd = [
        "git",
        "checkout",
        branch,
    ]
    subprocess.check_output(cmd, shell=False, stderr=subprocess.STDOUT)


def git_update_subtree(temp_repo_dir):
    """update any git subtree to the correct version.

    """
    print("Updating git subtrees....")
    externals_filename = "{0}/SVN_EXTERNAL_DIRECTORIES".format(temp_repo_dir)
    externals = []
    with open(externals_filename, 'r') as externals_file:
        externals = externals_file.readlines()

    for e in externals:
        ext = e.split()
        if len(ext) < 2:
            continue
        ext_dir = ext[0]
        url = ext[1]
        ext_commit = url.split('/')[-1]
        ext_url = '/'.join(url.split('/')[0:5])
        if ext_url.find('git') > 0:
            cmd = [
                'git',
                'subtree',
                'pull',
                '--prefix',
                ext_dir,
                ext_url,
                ext_commit,
            ]
            print("    {0}".format(' '.join(cmd)))
            subprocess.check_output(cmd, shell=False, stderr=subprocess.STDOUT)


def git_add_new_cesm(new_tag):
    """Add the new cesm files to git
    """
    print("Committing new cesm to git")
    cmd = [
        "git",
        "add",
        "-A",
    ]
    subprocess.check_output(cmd, shell=False, stderr=subprocess.STDOUT)

    cmd = [
        "git",
        "commit",
        "-m",
        "\'pull {0} tags from svn\'".format(new_tag)
    ]
    subprocess.check_output(cmd, shell=False, stderr=subprocess.STDOUT)


def git_status():
    """run the git status command
    """
    print("Running git status")
    cmd = [
        "git",
        "status",
    ]
    subprocess.check_output(cmd, shell=False, stderr=subprocess.STDOUT)


def push_to_origin_and_cleanup(branch, new_dir, temp_repo_dir):
    """
    """
    print("Pushing changes to git origin and removing update directory...")
    cmd = [
        "git",
        "push",
        "origin",
        branch,
    ]
    subprocess.check_output(cmd, shell=False, stderr=subprocess.STDOUT)
    os.chdir(new_dir)
    shutil.rmtree(temp_repo_dir)


# -------------------------------------------------------------------------------
#
# main
#
# -------------------------------------------------------------------------------

def main(options):

    config = read_config_file(options.config[0])
    new_tag = new_tag_from_config(config)

    # NOTE: just assume git is available in the path!
    cwd = os.getcwd()

    repo_dir = os.path.abspath("{0}/{1}".format(cwd, options.repo[0]))

    temp_repo_dir = "{0}/{1}-update-{2}".format(cwd, options.repo[0], new_tag)
    if os.path.isdir(temp_repo_dir):
        raise RuntimeError(
            "ERROR: temporary git repo dir already exists:\n    {0}".format(temp_repo_dir))

    clone_cesm_git(repo_dir, temp_repo_dir)
    os.chdir(temp_repo_dir)
    switch_git_branch(config["git"]["branch"])
    remove_current_working_copy()

    svn_checkout_cesm(config['cesm'], debug=options.debug)
    update_svn_externals(
        temp_repo_dir,
        config['cesm']['repo'],
        config['externals'])

    git_add_new_cesm(new_tag)
    git_update_subtree(temp_repo_dir)

    if (options.feelin_lucky):
        push_to_origin_and_cleanup(config["git"]["branch"], cwd, temp_repo_dir)

    print("Finished updating cesm to git.")
    return 0


if __name__ == "__main__":
    options = commandline_options()
    try:
        status = main(options)
        sys.exit(status)
    except Exception as error:
        print(str(error))
        if options.backtrace:
            traceback.print_exc()
        sys.exit(1)
