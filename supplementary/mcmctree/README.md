To reproduce, do the following:

```shell
$ cd run
$ mcmctree mcmctree.ctl
```

Then, after the first run finishes, there will be a file called `out.BV` inside 
the folder. Rename it to `in.BV` and edit the file `mcmctree.ctl`, changing the
option in line 7 to:

```
usedata = 2 in.BV
```

After that, save the configuration file and run MCMCTree again with the same
config file as input, inside the same directory.