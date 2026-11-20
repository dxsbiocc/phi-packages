# The spatial adjacency graph is the prerequisite for every spatial statistic:
# SVG, neighbourhood enrichment, co-occurrence and communication all read it.
# It is built on the tissue coordinates (obsm[coord_key]), NOT on the expression
# PCA/UMAP — two spots can be neighbours in tissue even when their expression
# differs. neighbour_method picks the geometry:
#   grid     — Visium hex/square lattice; n_neighs=6 matches the hexagonal layout
#   knn      — k nearest spots/bins; good for Stereo-seq bins or irregular layouts
#   radius   — everything within a physical distance; needs a meaningful radius
#   delaunay — triangulation; adaptive, no fixed k
# squidpy 1.8 deprecates the catch-all spatial_neighbors, so call the mode-specific
# builder directly. Results land in obsp['spatial_connectivities'/'spatial_distances'].
import squidpy as sq
from scipy.spatial import cKDTree

if coord_key not in adata.obsm:
    raise ValueError(f"coordinates '{coord_key}' not in obsm; spatial statistics need them")
if neighbor_method not in ("knn", "radius", "delaunay", "grid"):
    raise ValueError("neighbor_method must be knn, radius, delaunay or grid")

coords = np.asarray(adata.obsm[coord_key], dtype=float)
if neighbor_method == "grid":
    sq.gr.spatial_neighbors_grid(adata, spatial_key=coord_key, n_neighs=n_neighs, n_rings=1)
elif neighbor_method == "delaunay":
    sq.gr.spatial_neighbors_delaunay(adata, spatial_key=coord_key)
elif neighbor_method == "radius":
    r = radius
    if r is None:
        # Derive a radius from the typical spot spacing: 1.5x the median nearest
        # neighbour distance connects immediate neighbours without bridging gaps.
        dist, _ = cKDTree(coords).query(coords, k=2)
        r = float(np.median(dist[:, 1]) * 1.5)
        analysis_notes.append(f"radius derived from spacing: {r:.3g} coordinate units")
    sq.gr.spatial_neighbors_radius(adata, radius=r, spatial_key=coord_key)
else:
    sq.gr.spatial_neighbors_knn(adata, spatial_key=coord_key, n_neighs=n_neighs)

analysis_steps.append("spatial-neighbors")
graph = adata.obsp["spatial_connectivities"]
degree = np.asarray(graph.sum(axis=1)).ravel()
print(f"spatial graph ({neighbor_method}): {graph.nnz} edges, mean degree {degree.mean():.2f}")
if neighbor_method == "grid" and platform_info and platform_info.get("platform") != "visium":
    analysis_notes.append("grid adjacency assumes a regular lattice; for irregular bins prefer knn/radius")