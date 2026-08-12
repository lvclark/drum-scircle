make_umaps <- function(x, n = 10L, seeds = seq_len(n)){
    if(length(seeds) < n){
        stop("Need at least as many seeds as n.")
    }
    if(length(n) != 1){
        stop("Need just one value for n.")
    }
    out <- vector(mode = "list", length = n)

    for(i in seq_len(n)){
        set.seed(seeds[i])
        u <- umap::umap(x, method = "naive", preserve.seed = FALSE)
        out[[i]] <- u$layout
    }
    return(out)
}

rotate_umaps <- function(umap_list){
    out <- vector(mode = "list", length = length(umap_list))
    out[[1]] <- umap_list[[1]]
    for(i in seq_along(umap_list)[-1]){
        p <- pracma::procrustes(umap_list[[1]], umap_list[[i]])
        out[[i]] <- p$P
    }
    return(out)
}

order_umaps <- function(umap_list){
    n <- length(umap_list)
    # get distances between matrices
    d <- matrix(0, nrow = n, ncol = n)
    for(i in 1:(n - 1)){
        for(j in (i+1):n){
            d0 <- sum((umap_list[[i]] - umap_list[[j]]) ^ 2)
            d[i,j] <- d0
            d[j,i] <- d0
        }
    }
    # find a short path to loop through the umaps
    # don't exhaustively search for the shortest distance
    path <- integer(n)
    path[1] <- 1L
    last <- 1L
    i <- 2L
    while(any(path == 0L)){
        check <- seq_len(n)[-path[path != 0]]
        last <- check[which.min(d[last,check])]
        path[i] <- last
        i <- i + 1L
    }

    return(umap_list[path])
}

# Build long-format data frame for ggplot
prep_gg <- function(umap_list, color_by){
    n <- length(umap_list)
    df_list <- vector(mode = "list", length = n)
    for(i in seq_len(n)){
        df_list[[i]] <- data.frame(
            UMAP_1 = umap_list[[i]][,1],
            UMAP_2 = umap_list[[i]][,2],
            Color = color_by
        )
        df_list[[i]]$Iteration <- i
        if(is.null(rownames(umap_list[[i]]))){
            df_list[[i]]$Sample <-
                paste0("sample", seq_len(nrow(umap_list[[i]])))
        } else {
            df_list[[i]]$Sample <- rownames(umap_list[[i]])
        }
    }
    return(do.call(rbind, df_list))
}

# Build plot
plot_gg <- function(df){
    p <- ggplot2::ggplot(df, ggplot2::aes(x = UMAP_1, y = UMAP_2, color = Color)) +
        ggplot2::geom_point() +
        gganimate::transition_time(Iteration) +
        gganimate::ease_aes("linear")
    return(p)
}

# Testing
sample_meta <- read.csv("../../SIDS/ramirez_2025-12_downstream/popstruct/sample_metadata_2025-12-16.csv")
x <- as.matrix(sample_meta[,c("PC1", "PC2", "PC3")])

u <- make_umaps(x)
r <- rotate_umaps(u)
o <- order_umaps(r)
pr <- prep_gg(o, sample_meta$Q_African)

p <- plot_gg(pr) + labs(color = "African ancestry")

p

library(ggplot2)

ggplot(sample_meta, aes(x = u[[1]][,1], y = u[[1]][,2],
                        color = PC1)) +
    geom_point()

ggplot(sample_meta, aes(x = u[[2]][,1], y = u[[2]][,2],
                        color = PC1)) +
    geom_point()

ggplot(sample_meta, aes(x = u[[3]][,1], y = u[[3]][,2],
                        color = PC1)) +
    geom_point()

ggplot(sample_meta, aes(x = r[[3]][,1], y = r[[3]][,2],
                        color = PC1)) +
    geom_point()
