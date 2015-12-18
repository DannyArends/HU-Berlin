
graph <- rbind(
  c("A", "B"),
  c("B", "C"),
  c("B", "D"),
  c("E", "F"),
  c("F", "G"),
  c("A", "C")
)

getConnections <- function(graph, node){
  ci <- which(apply(graph, 1, function(x){any(x == node)} ))
  connectedTo <- unique(as.character(graph[ci, ]))
  connectedTo <- connectedTo[-which(connectedTo == node)]
  return(connectedTo)
}

traverseGraph <- function(graph, verbose = FALSE){
  nodes <- unique(as.character(graph))
  visited <- c()
  groups <- vector("list", length(nodes))
  group <- 1
  t <- 1
  while(length(nodes) > 0){
    cat("Group:", group)
    graph <- matrix(graph,,2)
    startNode <- names(which.max(table(as.character(graph))))
    cat(", start node", startNode, "Graph:", nrow(graph), "\n")
    connections <- getConnections(graph, startNode)
    groupmembers <- c(startNode)
    while(length(connections) > 0) {
      if(verbose) cat(length(connections), "\n")
      more <- getConnections(graph, connections[1])
      groupmembers <- unique(c(groupmembers, connections[1]))
      connections <- unique(c(connections[-1], more[which(!more %in% groupmembers)]))
    }
    cat("Groupsize:", length(groupmembers), "\n")
    toR <- apply(graph, 1, function(x){
      return(any(x %in% groupmembers))
    })
    groups[[group]] <- groupmembers
    graph <- graph[-which(toR),]
    nodes <- unique(as.character(graph))
    group <- group + 1
  }
  groups <- groups[1:(group-1)]
  return(groups)
}

traverseGraph(graph)