# Agreeable groups of higher genus

For a non-CM elliptic curve E over Q, the image of the adelic Galois representation of E is an open subgroup of GL(2,Zhat) whose agreeable closure we denote by G_E.

Assume that G_E is not contained in a group G with prime power level and X_G(Q) finite.

In a [previous project](https://github.com/davidzywina/OpenImage), up to conjugacy in GL(2,Zhat), we found that G_E must lie in one of finitely many agreeable subgroups G.  Many of the groups in the classification are only given implicitly.  The file `agreeable.dat` contains the information from this early work.

This code finds a complete sequence of such groups G.   We also look at some local obstructions to remove some groups G for which X_G(Q) is empty.  This list was made so that one could start determining the rational points of these modular curves X_G.

We make use of our faster modular curve code and the significantly faster group theory code of Drew Sutherland.

The list of groups can be found in the file `groups.m`.    The groups are produced by the field `find_groups.m` (it took a little over an hour on my machine; your times may vary).