# Track B prototype: spur-pruned, junction-merged branch tracing (the graph contract of
# Kovesi edgelink + cleanedgelist as used by MLV traceSkeleton, without its duplicated merge
# points and without dependence on edge enumeration order beyond documented tie keys).
trace_pruned <- function(skeleton, min_points = 10L) {
  h <- nrow(skeleton)
  paths <- ns$cpp_trace_graph(skeleton, 2L)
  E <- lapply(paths, as.integer); alive <- rep(TRUE, length(E))
  ends <- function(i) c(E[[i]][1], E[[i]][length(E[[i]])])
  degree <- function() { d <- table(unlist(lapply(which(alive), function(i) { e <- ends(i); if (e[1] == e[2] && length(E[[i]]) > 1) e else e })))
    setNames(as.integer(d), names(d)) }
  deg <- degree(); dg <- function(node) { v <- deg[as.character(node)]; ifelse(is.na(v), 0L, v) }
  len <- vapply(E, length, 1L)
  # 1. spurs: edges with one free end (degree 1) and the other at a junction (degree >= 3);
  #    at each junction remove only the shortest short spur (ties: smaller first pixel index).
  spur <- which(alive & len < min_points)
  info <- lapply(spur, function(i) { e <- ends(i); d <- dg(e); list(i = i, free = any(d == 1), junction = e[which.max(d)], jd = max(d), isolated = all(d == 1)) })
  remove <- integer()
  for (z in info) if (z$isolated) remove <- c(remove, z$i)
  cand <- Filter(function(z) z$free && !z$isolated && z$jd >= 3, info)
  if (length(cand)) {
    by_node <- split(cand, vapply(cand, function(z) as.character(z$junction), ""))
    for (g in by_node) { jd <- g[[1]]$jd
      if (jd >= 4) remove <- c(remove, vapply(g, `[[`, 1L, "i"))     # MLV: no merge needed, remove all
      else { o <- order(len[vapply(g, `[[`, 1L, "i")], vapply(g, function(z) E[[z$i]][1], 1L)); remove <- c(remove, g[[o[1]]]$i) } } }
  alive[unique(remove)] <- FALSE
  # 2. merge edges through nodes of degree exactly two (shared node kept once)
  repeat { deg <- degree(); two <- as.integer(names(deg)[deg == 2]); merged <- FALSE
    for (node in two) { inc <- which(alive & vapply(seq_along(E), function(i) alive[i] && any(ends(i) == node), TRUE))
      if (length(inc) != 2) next
      a <- E[[inc[1]]]; b <- E[[inc[2]]]
      if (a[length(a)] != node) a <- rev(a); if (b[1] != node) b <- rev(b)
      E[[inc[1]]] <- c(a, b[-1]); alive[inc[2]] <- FALSE; merged <- TRUE; break }
    if (!merged) break }
  keep <- which(alive); len <- vapply(E[keep], length, 1L)
  E <- E[keep]
  lapply(E, function(index) { rc <- arrayInd(index, dim(skeleton)); data.frame(row = rc[, 1], column = rc[, 2], x = rc[, 2], y = rc[, 1]) })
}
