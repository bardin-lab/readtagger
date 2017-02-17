Development
-----------

To run the test, install tox and then run tox.
This will execute pytest and planemo to handle all defined tests.

::

    pip install tox
    tox

Release
-------

Edit the HISTORY.rst file and bump the version:

::

    bumpversion minor --commit

This will change all version numbers.  Now push to master (or open a PR if you
don't have push rights). Do not push tags yet.  Make sure current master branch
builds correctly. Once all is good you can push the tags

::

    git push origin master --tags

If there is a problem in deployment, you can remove the current tag and retag a
different commit (here it's the latest commit):

::

    git push origin :refs/tags/<tag you want to remove>
    git tag -fa <tag you want to add>

