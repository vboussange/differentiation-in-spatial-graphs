All graphs are listed [here](https://users.cecs.anu.edu.au/~bdm/data/graphs.html)

To extract the graphs:

One needs to use the [gtools](https://users.cecs.anu.edu.au/~bdm/data/formats.html) package. 

```
./showg_mac64 graph7c.g6 -p1: -e -q -o1 -l0 > allgraphs_M7.txt
awk 'NR%2==0' allgraphs_M7.txt > allgraphs_M7_final.txt
sed -i .bak "s:  :),(:g" allgraphs_M7_final.txt
sed -i .bak -e 's/$/)]/' allgraphs_M7_final.txt
sed -i .bak -e 's/^/[(/' allgraphs_M7_final.txt
sed -i .bak "s: :,:g" allgraphs_M7_final.txt
rm allgraphs_M7.txt
```
