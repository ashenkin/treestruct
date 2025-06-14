# treestruct
R package for analysis and manipulation of tree structure models.  Wiki (generated by devin.ai) can be found here: https://deepwiki.com/ashenkin/treestruct/1-overview

This is an R package of tools for analyzing tree structural models based on terrestrial lidar data, and particularly those generated by Pasi Raumonen's treeqsm tool.

Note that the analyses of treestruct objects rely on the object attributes not being modified.  For example, it relies on internodes_reordered being correct about the states of the internodes in the object.  If these are set differently, results can be incorrect (e.g. pathlength calculations will be incorrect if internodes haven't actually been reordered).  The long and the short of it: only change the attributes directly if you know what you're doing.

Collaborators are encouraged!  If you'd like to contribute, please fork the repo and submit pull requests with changes.

Example usage:

```r
library(treestruct)
my_trees = import_treestructs_from_dir(qsm_path = "/path/to/cyl/files/"), qsmver = "UCL", recursive = T)
my_qsms = TreeStructs(dataset = "my_trees", treestructs = my_trees)
my_qsms = run_all(my_qsms)  # run all architectural analyses
my_qsms = parse_id(my_qsms, treetag_regex = ".*(?=\\.mat)", nobranchcode = T)  # set tree tag in the data structure
my_treestructs = getTreestructs(my_qsms)  # take a look at the resulting data
```

This software is released under the CC-BY-NC-SA-4.0 license.  See https://creativecommons.org/licenses/by-nc-sa/4.0/ .
