using Graphs, GLMakie, GraphMakie
import GraphMakie: Spring
GLMakie.activate!()

pts = rand(3, 10) # 10 vertices in R^3

g, dists = euclidean_graph(pts, p=1, bc=:periodic) # Taxicab-distance (L^1);

e = edges(g)

graphplot(g; layout=GraphMakie.Spring(dim=3))


