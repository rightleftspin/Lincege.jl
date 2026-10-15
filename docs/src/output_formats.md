# Output Formats

`write_to_json` and `write_to_txt_file` export the same cluster data in two
formats. 

## Single expansion

Use either writer with the same expansion and lattice:

```julia
write_to_json(expansion, lattice, "clusters.json")
write_to_txt_file(expansion, lattice, "clusters.txt")
```

### JSON

The file is an array of cluster objects:

```json
[
  {
    "cluster_hash": "###",
    "order": 3,
    "n_sites": 3,
    "coordinates": [[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]],
    "site_colors": [1, 1, 1],
    "bonds": [[1, 2, 1], [2, 3, 1]],
    "weights": [0.0, 0.0, 6.0]
  }
]
```

### Text

The first line names the fields. Each subsequent line represents one cluster:

```text
# n_sites : cluster_hash : bonds : coordinates : site_colors : weights : order
3 : 100 : 1 2 1 2 3 1 : (0.0,0.0) (1.0,0.0) (2.0,0.0) : 1 1 1 : 0.0 0.0 6.0 : 3
```

## Joint isomorphic and symmetric expansions

Pass both expansions and their cluster sets, built on the same lattice:

```julia
write_to_json(
        isomorphic_expansion, isomorphic_clusters,
        symmetric_expansion, symmetric_clusters,
        lattice, "related_clusters.json",
)
write_to_txt_file(
        isomorphic_expansion, isomorphic_clusters,
        symmetric_expansion, symmetric_clusters,
        lattice, "related_clusters.txt",
)
```

### JSON

Each isomorphic cluster contains a `related_symmetric` array. In this example,
a straight chain and a bent chain share the same graph connectivity:

```json
[
  {
    "cluster_hash": "###",
    "order": 3,
    "n_sites": 3,
    "canonical_bonds": [[1, 2, 1], [2, 3, 1]],
    "canonical_site_colors": [1, 1, 1],
    "weights": [0.0, 0.0, 6.0],
    "related_symmetric": [
      {
        "cluster_hash": "####",
        "permuted_coordinates": [[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]],
        "weights": [0.0, 0.0, 2.0]
      },
      {
        "cluster_hash": "#####",
        "permuted_coordinates": [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]],
        "weights": [0.0, 0.0, 4.0]
      }
    ]
  }
]
```

### Text

The joint text file has two header lines and two record types:

```text
# isomorphic : n_sites : cluster_hash : canonical_bonds : canonical_site_colors : weights : order
# symmetric : isomorphic_hash : cluster_hash : permuted_coordinates : weights
isomorphic : 3 : 100 : 1 2 1 2 3 1 : 1 1 1 : 0.0 0.0 6.0 : 3
symmetric : 100 : 201 : (0.0,0.0) (1.0,0.0) (2.0,0.0) : 0.0 0.0 2.0
symmetric : 100 : 202 : (0.0,0.0) (1.0,0.0) (1.0,1.0) : 0.0 0.0 4.0
```
