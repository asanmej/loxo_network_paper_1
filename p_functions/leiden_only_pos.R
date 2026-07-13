# Author: Alejandro Santos Mejías
# Date last update: 2026-06-12
# Script name: leiden_only_pos.R
# Aims: To run Leiden algorithm

leiden_only_pos <- function(edges, nodes, nombre_base, ruta_subcarpeta_estrato){
  
  if (!dir.exists(ruta_subcarpeta_estrato)) {
    dir.create(ruta_subcarpeta_estrato, recursive = TRUE)
  }
  
  grafo_X <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)
  
  # ============================================================================
  # 5. OPTIMIZACIÓN DE RESOLUCIÓN (LEIDEN)
  # ============================================================================
  res_values <- seq(0.005, 1.5, by = 0.001)
  resultados <- data.frame(resolucion = numeric(), modularidad = numeric(), num_clusters = integer())
  
  for (res in res_values) {
    # IMPORTANTE: el parámetro correcto en igraph es 'resolution_parameter'
    clustering <- cluster_leiden(grafo_X, objective_function = "modularity", resolution_parameter = res)
    # He tenido que volver a poner los pesos en absoluto la red los ha vuelto a transformar en negativo 
    mod <- modularity(grafo_X, membership = clustering$membership, weights = E(grafo_X)$weight)
    
    num_clusters <- length(unique(clustering$membership))
    
    resultados <- resultados %>%
      add_row(resolucion = res, modularidad = mod, num_clusters = num_clusters)
  }
  
  # Extraer la mejor resolución
  mejor <- resultados %>%
    filter(modularidad == max(modularidad)) %>%
    slice(1)
  
  cat("\n🟢 Mejor resolución:", mejor$resolucion, 
      "con modularidad =", round(mejor$modularidad, 4),
      "y número de clústeres =", mejor$num_clusters, "\n")
  
  # ============================================================================
  # 6 y 7. APLICAR COMUNIDADES AL GRAFO FINAL
  # ============================================================================
  clustering_final <- cluster_leiden(grafo_X, objective_function = "modularity", resolution_parameter = mejor$resolucion, n_iterations = 100)
  
  mejor$num_clusters == max(clustering_final$membership)
  
  nodes$community <- factor(clustering_final$membership)
  V(grafo_X)$size <- nodes$size
  V(grafo_X)$community <- factor(clustering_final$membership)
  
  membership_vector <- as.numeric(V(grafo_X)$community)
  mod_final <- mod <- modularity(grafo_X, membership = clustering_final$membership, weights = E(grafo_X)$weight)
  cat("La modularidad es:", mod_final, "\n")
  
  # ============================================================================
  # 8 al 12. ESTÉTICA Y VISUALIZACIÓN
  # ============================================================================
  # Asignación de colores
  paleta_nivel <- rainbow(length(unique(membership_vector)))
  colores_asignados <- paleta_nivel[as.numeric(factor(membership_vector))]
  
  # Actualización del dataframe de nodos incluyendo la comunidad para visOptions
  nodes_nivel <- nodes %>%
    mutate(
      color = colores_asignados,
      community = membership_vector,
      group = factor(membership_vector),
      title = paste0(
        "<B>",
        label,
        "</B><br>Prevalencia: ",
        round(prevalencia_pct, 2),
        "%<br>Comunidad: ",
        membership_vector
      )
    )
  
  # Creamos el objeto igraph para este nivel
  grafo_nivel <- graph_from_data_frame(d = edges ,
                                       vertices = nodes_nivel,
                                       directed = FALSE)
  
  # Asignación de atributos de aristas
  E(grafo_nivel)$width <- abs(E(grafo_nivel)$weight) * 10
  E(grafo_nivel)$color <- ifelse(E(grafo_nivel)$correlation >= 0, "#1f78b4", "#e31a1c")
  
  # Generación del gráfico usando visIgraph
  vis_plot_nivel <- visIgraph(grafo_nivel, layout = "layout_with_fr") %>%
    visOptions(
      highlightNearest = list(enabled = TRUE, hover = TRUE),
      nodesIdSelection = TRUE,
      selectedBy = "community"
    ) %>%
    visPhysics(
      enabled = TRUE,
      solver = "forceAtlas2Based",
      forceAtlas2Based = list(
        gravitationalConstant = -50,
        springLength = 100
      )
    ) %>%
    htmlwidgets::prependContent(
      htmltools::tags$style(
        "div.visNetwork { width: 100% !important; height: 92vh !important; }"
      )
    )
  
  # Guardado automático en la subcarpeta
  nombre_archivo <- paste0("grafo_Leiden_", nombre_base, ".html")
  htmlwidgets::saveWidget(
    vis_plot_nivel,
    file.path(ruta_subcarpeta_estrato, nombre_archivo),
    selfcontained = TRUE
  )
  
  resumen_sbm <- nodes_nivel %>%
    group_by(group) %>%
    summarize(
      num_elementos = n(),
      prevalencia_media_bloque = mean(prevalencia_pct),
      enfermedades = paste(label, collapse = " | ")
    )
  write_csv(resumen_sbm, file.path(ruta_subcarpeta_estrato, paste0("resumen_bloques_Leiden_", nombre_base, ".csv")))
  
  # ============================================================================
  # 13. GUARDAR LOS RESULTADOS
  # ============================================================================
  sufix <- "_ebicglasso0.25"
  # Guardamos los nodos en la subcarpeta del estrato correspondiente
  write_csv(nodes, file.path(ruta_subcarpeta_estrato, paste0("nodos_comunidades_", nombre_base, sufix, ".csv")))
  
  # Guardamos la matriz de aristas filtrada en la subcarpeta del estrato correspondiente
  write_csv(edges, file.path(ruta_subcarpeta_estrato, paste0("aristas_filtradas_", nombre_base, sufix, ".csv")))
  
  cat("💾 Resultados de esta iteración guardados con éxito en:", ruta_subcarpeta_estrato, "\n")
  
  
}