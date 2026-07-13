# Author: Alejandro Santos Mejías
# Date last update: 2026-06-12
# Script name: spinnglass_pos.R
# Aims: To run Spinnglass algorithm without negative correlations.

spinnglass_pos <- function(edges, nodes, nombre_base, ruta_subcarpeta_estrato){
  
  if (!dir.exists(ruta_subcarpeta_estrato)) {
    dir.create(ruta_subcarpeta_estrato, recursive = TRUE)
  }
  
  grafo_X <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)
  
  # ============================================================================
  # 4. OPTIMIZACIÓN DE PARÁMETROS (SPINGLASS)
  # ============================================================================
  # El algoritmo Spin Glass requiere una estructura conectada. Si existen subcomponentes 
  # aisladas debido al umbral estricto, el bucle se ejecuta sobre la componente principal.
  if (!is_connected(grafo_X)) {
    componentes <- components(grafo_X)
    mayor_componente <- which.max(componentes$csize)
    nodos_filtrados <- V(grafo_X)$name[componentes$membership == mayor_componente]
    grafo_optimizacion <- induced_subgraph(grafo_X, nodos_filtrados)
  } else {
    grafo_optimizacion <- grafo_X
  }
  
  # Exploración del parámetro gamma característico del modelo de espines
  gamma_values <- seq(0.5, 2.0, by = 0.1)
  resultados <- tibble(gamma_spinglass = numeric(), modularidad = numeric(), num_clusters = integer())
  
  for (g in gamma_values) {
    tryCatch({
      clustering <- cluster_spinglass(grafo_optimizacion, weights = E(grafo_optimizacion)$weight, gamma = g)
      mod <- modularity(grafo_optimizacion, membership = clustering$membership, weights = E(grafo_optimizacion)$weight)
      num_clusters <- length(unique(clustering$membership))
      
      resultados <- resultados %>%
        add_row(gamma_spinglass = g, modularidad = mod, num_clusters = num_clusters)
    }, error = function(e) {
      # Control de excepciones por si alguna configuración de gamma genera inestabilidad numérica
    })
  }
  
  # Selección del parámetro gamma que maximiza la modularidad del sistema
  mejor <- resultados %>%
    filter(modularidad == max(modularidad)) %>%
    slice(1)
  
  cat("\n🟢 Mejor gamma Spinglass:", mejor$gamma_spinglass, 
      "con modularidad =", round(mejor$modularidad, 4),
      "y número de clústeres =", mejor$num_clusters, "\n")
  
  # ============================================================================
  # 5. APLICAR COMUNIDADES AL GRAFO FINAL
  # ============================================================================
  clustering_final <- cluster_spinglass(grafo_optimizacion, weights = E(grafo_optimizacion)$weight, gamma = mejor$gamma_spinglass)
  
  # Mapeo de las pertenencias obtenidas al listado general de nodos
  nodos_con_comunidad <- tibble(
    id = V(grafo_optimizacion)$name,
    community = as.numeric(clustering_final$membership)
  )
  
  nodes <- nodes %>%
    left_join(nodos_con_comunidad, by = "id") %>%
    mutate(community = ifelse(is.na(community), 0, community)) 
  
  nodes$community <- factor(nodes$community)
  V(grafo_X)$size <- nodes$size
  V(grafo_X)$community <- nodes$community
  
  membership_vector <- as.numeric(V(grafo_X)$community)
  
  # Recalcular la modularidad del grafo optimizado para la salida en consola
  mod_final <- modularity(grafo_optimizacion, membership = clustering_final$membership, weights = E(grafo_optimizacion)$weight)
  cat("La modularidad es:", mod_final, "\n")
  
  # ============================================================================
  # 6. ESTÉTICA Y VISUALIZACIÓN INTERACTIVA
  # ============================================================================
  paleta_nivel <- rainbow(length(unique(membership_vector)))
  colores_asignados <- paleta_nivel[as.numeric(factor(membership_vector))]
  
  nodes_nivel <- nodes %>%
    mutate(
      color = colores_asignados,
      community = as.numeric(community),
      group = factor(community),
      title = paste0(
        "<B>",
        label,
        "</B><br>Prevalencia: ",
        round(prevalencia_pct, 2),
        "%<br>Comunidad: ",
        community
      )
    )
  
  grafo_nivel <- graph_from_data_frame(d = edges, vertices = nodes_nivel, directed = FALSE)
  
  # Configuración visual de las conexiones según peso y signo
  E(grafo_nivel)$width <- abs(E(grafo_nivel)$weight) * 10
  E(grafo_nivel)$color <- ifelse(E(grafo_nivel)$correlation >= 0, "#1f78b4", "#e31a1c")
  
  # Renderizado de la red con el motor visNetwork utilizando ForceAtlas2
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
  
  # Guardar el archivo HTML resultante en la ruta del estrato
  nombre_archivo <- paste0("grafo_Spinglass_", nombre_base, ".html")
  htmlwidgets::saveWidget(
    vis_plot_nivel,
    file.path(ruta_subcarpeta_estrato, nombre_archivo),
    selfcontained = TRUE
  )
  
  # Crear y exportar el informe resumido por cada clúster/bloque identificado
  resumen_sbm <- nodes_nivel %>%
    group_by(group) %>%
    summarize(
      num_elementos = n(),
      prevalencia_media_bloque = mean(prevalencia_pct),
      enfermedades = paste(label, collapse = " | ")
    )
  write_csv(resumen_sbm, file.path(ruta_subcarpeta_estrato, paste0("resumen_bloques_Spinglass_", nombre_base, ".csv")))
  
  # ============================================================================
  # 7. EXPORTACIÓN DE TABLAS DE RESULTADOS
  # ============================================================================
  sufix <- "_ebicglasso0.25"
  write_csv(nodes, file.path(ruta_subcarpeta_estrato, paste0("nodos_comunidades_", nombre_base, sufix, ".csv")))
  write_csv(edges, file.path(ruta_subcarpeta_estrato, paste0("aristas_filtradas_", nombre_base, sufix, ".csv")))
  
  cat("💾 Resultados de esta iteración guardados con éxito en:", ruta_subcarpeta_estrato, "\n")
  
  
}