# Author: Alejandro Santos Mejías
# Date last update: 2026-06-12
# Script name: spinglass_neg.R
# Aims: To run negative implementation of Spinnglass algorithm


spinnglass_neg <- function(edges, nodes, nombre_base, ruta_subcarpeta_estrato) {
  
  if (!dir.exists(ruta_subcarpeta_estrato)) {
    dir.create(ruta_subcarpeta_estrato, recursive = TRUE)
  }
  
  grafo_X <- graph_from_data_frame(d = edges,
                                   vertices = nodes,
                                   directed = FALSE)
  
  # ============================================================================
  # 4. OPTIMIZACIÓN DE PARÁMETROS (SPINGLASS CON SIGNOS)
  # ============================================================================
  # Spinglass con enlaces negativos requiere inspeccionar la conectividad basal.
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
  n_replicas <- 50 # Número de simulaciones por gamma para evaluar estabilidad
  resultados_estabilidad <- tibble(
    gamma_spinglass = numeric(),
    ami_promedio = numeric(),
    num_clusters_mediana = numeric()
  )
  
  for (g in gamma_values) {
    # Matriz para guardar las membresías de cada réplica de este nivel de gamma
    # Filas: nodos del grafo, Columnas: réplicas
    membresias_gamma <- matrix(NA, nrow = vcount(grafo_optimizacion), ncol = n_replicas)
    clusters_detectados <- c()
    
    for (r in 1:n_replicas) {
      tryCatch({
        clustering <- cluster_spinglass(
          grafo_optimizacion,
          weights = E(grafo_optimizacion)$correlation,
          gamma = g,
          implementation = "neg"
        )
        membresias_gamma[, r] <- clustering$membership
        clusters_detectados <- c(clusters_detectados, length(unique(clustering$membership)))
      }, error = function(e) {
        # Control de excepciones por si el recocido simulado diverge numéricamente
        message(paste0('ERROR: Divergencia en gamma = ', g))
      })
    }
    
    # Eliminamos columnas vacías si alguna réplica falló por error numérico
    membresias_gamma <- membresias_gamma[, colSums(is.na(membresias_gamma)) == 0, drop = FALSE]
    
    # Si tenemos suficientes réplicas válidas, calculamos la consistencia interna
    if (ncol(membresias_gamma) > 1) {
      valores_ami <- c()
      
      # Comparamos todos los pares únicos de particiones para medir la estabilidad
      for (i in 1:(ncol(membresias_gamma) - 1)) {
        for (j in (i + 1):ncol(membresias_gamma)) {
          # Uso exclusivo de igraph::compare con el método de Información Mutua Ajustada
          ami_par <- compare(membresias_gamma[, i], membresias_gamma[, j], method = "nmi")
          valores_ami <- c(valores_ami, ami_par)
        }
      }
      
      resultados_estabilidad <- resultados_estabilidad %>%
        add_row(
          gamma_spinglass = g,
          ami_promedio = mean(valores_ami, na.rm = TRUE),
          num_clusters_mediana = median(clusters_detectados, na.rm = TRUE)
        )
    }
  }
  
  # Selección del parámetro gamma que maximiza la estabilidad del sistema (AMI promedio)
  mejor_estabilidad <- resultados_estabilidad %>%
    filter(ami_promedio == max(ami_promedio, na.rm = TRUE)) %>%
    slice(1)
  
  cat(
    "\n🟢 Mejor gamma Spinglass por Criterio de Estabilidad (AMI):",
    mejor_estabilidad$gamma_spinglass,
    "\nEstabilidad (AMI Promedio) =",
    mejor_estabilidad$ami_promedio,
    "\nNúmero mediano de clústeres =",
    mejor_estabilidad$num_clusters_mediana,
    "\n"
  )
  
  # ============================================================================
  # 5. APLICAR COMUNIDADES AL GRAFO FINAL
  # ============================================================================
  # Pasamos el vector con signos correspondientes a 'correlation' para la segmentación definitiva
  clustering_final <- cluster_spinglass(
    grafo_optimizacion,
    weights = E(grafo_optimizacion)$correlation,
    gamma = mejor_estabilidad$gamma_spinglass,
    implementation = "neg"
  )
  
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
  
  # Recalcular la modularidad usando el vector con signo original
  # mod_final <- modularity(grafo_optimizacion, membership = clustering_final$membership, weights = E(grafo_optimizacion)$correlation)
  # cat("La modularidad es:", mod_final, "\n")
  
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
  
  grafo_nivel <- graph_from_data_frame(d = edges,
                                       vertices = nodes_nivel,
                                       directed = FALSE)
  
  # Configuración visual de las conexiones utilizando el valor absoluto para el grosor
  E(grafo_nivel)$width <- E(grafo_nivel)$weight * 10
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
  write_csv(resumen_sbm, file.path(
    ruta_subcarpeta_estrato,
    paste0("resumen_bloques_Spinglass_", nombre_base, ".csv")
  ))
  
  # ============================================================================
  # 7. EXPORTACIÓN DE TABLAS DE RESULTADOS
  # ============================================================================
  sufix <- "_ebicglasso0.25"
  write_csv(nodes, file.path(
    ruta_subcarpeta_estrato,
    paste0("nodos_comunidades_", nombre_base, sufix, ".csv")
  ))
  write_csv(edges, file.path(
    ruta_subcarpeta_estrato,
    paste0("aristas_filtradas_", nombre_base, sufix, ".csv")
  ))
  
  cat(
    "💾 Resultados de esta iteración guardados con éxito en:",
    ruta_subcarpeta_estrato,
    "\n"
  )
}