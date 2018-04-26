
library("igraph")
library("network")
library("sna")
library("ndtv")

nodes<-read.csv("./networks/hawaii_vertices.csv")
stepedges<-read.csv("./networks/step_edges.csv",header=F)
colrs<-c("black","gray50")

#island model make all possible connections

edges<-t(combn(x=nodes$vertex,m=2))
# and the reverse
edges2<-rbind(edges,edges[,c(2,1)])

#make the n-island model
island<-graph.data.frame(edges2,nodes,directed=T)

circle<-layout_in_circle(island,order=c(11,10,9,8,7,6,5,4,3,2,1,14,13,12))
hi_islands<-as.matrix(nodes[,c(7,6)])

plot(island,layout=hi_islands, edge.arrow.size=0.05, edge.curved=0.3, vertex.label=V(island)$label, vertex.color=V(island)$sector, vertex.size=8,vertex.label.cex=0.7,vertex.shape="csquare")


#make a stepping-stone model
stepstone<-graph.data.frame(stepedges,nodes,directed=T)

plot(stepstone,layout=hi_islands, edge.arrow.size=0.05, edge.curved=0.3, vertex.label=V(stepstone)$label, vertex.color=V(stepstone)$sector, edge.color="black",vertex.label.cex=0.7,vertex.size=8)

plot(stepstone,layout=hi_islands, edge.arrow.size=0.05, edge.curved=0.3, vertex.label=V(stepstone)$label, vertex.color=V(stepstone)$sector, edge.color="gray",vertex.label.cex=0.7,vertex.size=8, vertex.shape="csquare")


# make simulated graphs
stepedges<-stepedges[1:18,]
nodes<-nodes[1:10,]

edges<-t(combn(x=nodes$vertex,m=2))
# and the reverse
edges2<-rbind(edges,edges[,c(2,1)])

island<-graph.data.frame(edges2,nodes,directed=T)


oneD_lattice<-cbind(x=rep(1,10), y=seq(1,10,1))

island<-graph.data.frame(edges2,nodes,directed=T)

plot(island,layout=oneD_lattice, edge.arrow.size=0.05, edge.curved=0.3, vertex.color="grey", vertex.label=NA,vertex.size=8,vertex.label.cex=0.7,vertex.shape="csquare")


stepstone<-graph.data.frame(stepedges,nodes,directed=T)

plot(stepstone,layout=oneD_lattice, edge.arrow.size=0.05, edge.curved=0.3, vertex.label=NA, vertex.color="grey", edge.color="black",vertex.label.cex=0.7,vertex.size=8)

plot(stepstone,layout=oneD_lattice, edge.arrow.size=0.05, edge.curved=0.3, vertex.label=NA, vertex.color="grey", edge.color="grey",vertex.label.cex=0.7,vertex.size=8, vertex.shape="csquare")
