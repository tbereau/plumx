# Pitfall: Check H-bond interaction #

The major pitfall in using the code is to fail/forget to activate the H-bond interaction. One can easily check whether it is properly running by simulating the decaalanine example:
```
/PATH/TO/plumx/decaalanine
```

Follow the instructions of [the peptide tutorial](helloWorldPep.md) applied to the decaalanine example.

  1. Run it first **without** setting the GMX\_NBLISTCG environment variable. A 20-ps-long simulation (default) should be enough. The helix will quickly unfold.
  1. Set and export GMX\_NBLISTCG (following the [instructions](helloWorldPep.md)). The helix will be stable.