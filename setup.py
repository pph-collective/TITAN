import setuptools
import subprocess

# ========= FROM https://github.com/vivin/better-setuptools-git-version/blob/master/better_setuptools_git_version.py with edits =========
def get_latest_tag():
    return subprocess.getoutput("git describe --tags `git rev-list --tags --max-count=1`")

def get_tag_commit_sha(tag):
    """Return the commit that the tag is pointing to."""
    return subprocess.getoutput("git rev-list -n 1 {tag}".format(tag=tag))

def get_head_sha():
    """Return the sha key of HEAD."""
    return subprocess.getoutput('git rev-parse HEAD')

def is_head_at_tag(tag):
    """Return True or False depending on whether the given tag is pointing to HEAD"""
    return get_head_sha() == get_tag_commit_sha(tag)

def get_version():
    """
    Return the full git version using the given template. If there are no annotated tags, the version specified by
    starting_version will be used. If HEAD is at the tag, the version will be the tag itself. If there are commits ahead
    of the tag, the first 8 characters of the sha of the HEAD commit will be included.
    In all of the above cases, if the working tree is also dirty or contains untracked files, a "+dirty" suffix will be
    appended to the version.
    Returns:
        the formatted version based on tags in the git repository.
    """
    template="{tag}.dev{sha}"
    tag = get_latest_tag()
    if is_head_at_tag(tag):
        version = tag
    else:
        sha = get_head_sha()[:8]
        version = template.format(tag=tag, sha=sha)

    return version

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open('requirements.txt') as f:
    required = f.read().splitlines()

setuptools.setup(
    name="titan", # Replace with your own username,
    version=get_version(),
    # author="Example Author",
    # author_email="author@example.com",
    description="TITAN Agent Based Model",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/marshall-lab/TITAN",
    packages=setuptools.find_packages(),
    install_requires=required,
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    package_data={'titan': ['params/*.yml', 'settings/*/*.yml']},
    entry_points={
        'console_scripts': ['run_titan=titan:script_init']
    }
)
