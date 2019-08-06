# Announcement: The contrib-snufflings repository is leaving GitHub

*Potsdam, 2019-08-05*

Since last week, [GitHub is restricting access to their services based on
user nationality and residence](https://help.github.com/en/articles/github-and-trade-controls>) ([see
also](https://techcrunch.com/2019/07/29/github-ban-sanctioned-countries)).
Such restrictions are incompatible with scientific standards in
international research communities like seismology.

The contrib-snufflings software package is used by researchers worldwide. As researchers, we are obligated to retain open
access to all. To achieve this, we are now migrating our code repositories
away from GitHub to a new safe home. The new home of the contrib-snufflings repository
is at [git.pyrocko.org](https://git.pyrocko.org/pyrocko/contrib-snufflings/), open now.

To ensure a smooth
transition, we will keep a read-only version of the contrib-snufflings repository
at GitHub until 2019-10-01, when it will be deleted.

To update the upstream url of a cloned contrib-snufflings repository, run

```
git remote set-url origin https://git.pyrocko.org/pyrocko/contrib-snufflings.git
```

in the cloned directory.

To obtain a fresh clone, run

```
git clone https://git.pyrocko.org/pyrocko/contrib-snufflings.git contrib-snufflings
```

Thanks to the worldwide seismology community for all the support and help.

Best regards

*The contrib-snufflings Developers*
