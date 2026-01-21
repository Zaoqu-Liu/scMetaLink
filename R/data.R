#' @title CRC Example Expression Data
#' @description Example single-cell RNA-seq expression matrix from colorectal cancer
#'
#' @format A sparse dgCMatrix with 4,210 genes (rows) x 2,850 cells (columns).
#'   Only genes present in MetalinksDB are included to reduce file size.
#'
#' @details
#' This is a subset of colorectal cancer single-cell data containing:
#' \itemize{
#'   \item Tumor cells: Tumor Epithelial (600 cells)
#'   \item Immune cells: T (500), Plasma (250), B (150), TAM (150), Monocyte (120),
#'         Normal Macrophage (150), Mast (30)
#'   \item Stromal cells: CAF (200), Normal Fibroblast (200), Endothelial (100),
#'         Pericyte (50), SMC (30)
#'   \item Epithelial: Normal Epithelial (300)
#'   \item Other: Gliacyte (20)
#' }
#'
#' @source CellScope package example data
#'
#' @examples
#' \donttest{
#' data(crc_expr)
#' data(crc_meta)
#'
#' # Run scMetaLink analysis
#' obj <- createScMetaLink(crc_expr, crc_meta, "cell_type")
#' }
"crc_expr"

#' @title CRC Example Cell Metadata
#' @description Cell metadata for the CRC example dataset
#'
#' @format A data.frame with 2,850 rows (cells) and 3 columns:
#' \describe{
#'   \item{cell_type}{Cell type annotation (15 types)}
#'   \item{tumor_normal}{Tumor or Normal tissue origin}
#'   \item{tissue_region}{Tissue region}
#' }
#'
#' @examples
#' \donttest{
#' data(crc_meta)
#' table(crc_meta$cell_type)
#' }
"crc_meta"

#' @title MetalinksDB Database
#' @description Pre-compiled metabolite-protein interaction database from MetalinksDB
#'
#' @format A list containing the following components:
#' \describe{
#'   \item{edges}{data.frame with 41,894 metabolite-protein interactions containing:
#'     \itemize{
#'       \item hmdb: HMDB metabolite identifier
#'       \item uniprot: UniProt protein identifier
#'       \item source: Data source
#'       \item combined_score: Interaction confidence score (0-1000)
#'       \item mor: Mode of regulation (1=producing, -1=degrading, 0=binding)
#'       \item type: Interaction type ("lr"=ligand-receptor, "pd"=produce-degrade)
#'       \item transport_direction: For transporters ("in" or "out")
#'     }
#'   }
#'   \item{metabolites}{data.frame with 1,128 metabolites containing:
#'     \itemize{
#'       \item hmdb: HMDB identifier
#'       \item metabolite: Metabolite name
#'       \item pubchem: PubChem identifier
#'       \item metabolite_subclass: Chemical classification
#'     }
#'   }
#'   \item{proteins}{data.frame with 4,374 proteins containing:
#'     \itemize{
#'       \item uniprot: UniProt identifier
#'       \item gene_symbol: Gene symbol
#'       \item protein_type: Protein classification (enzyme, gpcr, transporter, etc.)
#'     }
#'   }
#'   \item{pathway}{data.frame with 157,741 metabolite-pathway associations}
#'   \item{cell_location}{data.frame with 2,816 subcellular location annotations}
#'   \item{tissue_location}{data.frame with 2,410 tissue location annotations}
#'   \item{disease}{data.frame with 3,216 disease associations}
#' }
#'
#' @details
#' The database enables two types of metabolite-mediated communication inference:
#'
#' **Ligand-Receptor (lr) type**: Direct metabolite-receptor binding interactions,
#' primarily involving GPCRs, nuclear hormone receptors, and ion channels.
#'
#' **Produce-Degrade (pd) type**: Enzyme-mediated metabolite production and
#' consumption, enabling inference of metabolite availability from enzyme expression.
#'
#' @source MetalinksDB (https://metalinks.org/)
#' @references
#' Schafer, S., et al. (2023). MetalinksDB: a knowledgebase of metabolite-centric
#' signaling. Nature Communications.
#'
#' @examples
#' \donttest{
#' # Access the database
#' db <- scMetaLink:::.load_metalinksdb()
#'
#' # View available metabolites
#' head(db$metabolites)
#'
#' # Check interaction types
#' table(db$edges$type)
#' }
"metalinksdb"

#' Load MetalinksDB Database
#'
#' @return A list containing the MetalinksDB data
#' @keywords internal
.load_metalinksdb <- function() {
  # Check if already loaded in package environment
  pkg_env <- getNamespace("scMetaLink")

  if (exists("metalinksdb", envir = pkg_env)) {
    return(get("metalinksdb", envir = pkg_env))
  }

  # Load from internal data
  data("metalinksdb", package = "scMetaLink", envir = environment())
  metalinksdb
}

#' Clean protein_type field (remove extra quotes)
#' @keywords internal
.clean_protein_type <- function(type_str) {
  if (is.na(type_str) || is.null(type_str)) {
    return(NA_character_)
  }
  # Remove surrounding quotes if present
  gsub('^"|"$', "", type_str)
}

#' Get Production Enzymes for Metabolites
#'
#' @description Extract enzymes and secretion transporters involved in metabolite
#'   production and release. This function filters for proteins that contribute to
#'   making metabolites available for intercellular communication.
#'
#' @details
#' The production potential consists of two components:
#' \itemize{
#'   \item Synthesis enzymes: Proteins that catalyze metabolite biosynthesis
#'   \item Secretion transporters: Membrane transporters that release metabolites
#'     to the extracellular space (transport_direction = "out")
#' }
#'
#' Uptake transporters (transport_direction = "in") are excluded from production
#' as they represent metabolite sensing/uptake rather than production/release.
#'
#' @param db MetalinksDB database list
#' @return data.frame of production relationships
#' @keywords internal
.get_production_enzymes <- function(db) {
  # Filter for production relationships (mor = 1, type = pd)
  prod_edges <- db$edges[db$edges$type == "pd" & db$edges$mor == 1, ]


  # Exclude uptake transporters (direction = "in") from production

  # Rationale: Uptake represents sensing/receiving, not production/release

  # Keep: enzymes (NA direction) and secretion transporters (direction = "out")
  prod_edges <- prod_edges[
    is.na(prod_edges$transport_direction) |
      prod_edges$transport_direction == "out",
  ]

  # Merge with gene symbols
  prod_edges <- merge(prod_edges, db$proteins[, c("uniprot", "gene_symbol")],
    by = "uniprot", all.x = TRUE
  )

  # Merge with metabolite names
  prod_edges <- merge(prod_edges, db$metabolites[, c("hmdb", "metabolite")],
    by = "hmdb", all.x = TRUE
  )

  prod_edges
}

#' Get Degradation Enzymes for Metabolites
#'
#' @param db MetalinksDB database list
#' @return data.frame of degradation relationships
#' @keywords internal
.get_degradation_enzymes <- function(db) {
  # Filter for degradation relationships (mor = -1, type = pd)
  deg_edges <- db$edges[db$edges$type == "pd" & db$edges$mor == -1, ]

  # Merge with gene symbols
  deg_edges <- merge(deg_edges, db$proteins[, c("uniprot", "gene_symbol")],
    by = "uniprot", all.x = TRUE
  )

  # Merge with metabolite names
  deg_edges <- merge(deg_edges, db$metabolites[, c("hmdb", "metabolite")],
    by = "hmdb", all.x = TRUE
  )

  deg_edges
}

#' Get Receptor/Sensing Proteins for Metabolites
#'
#' @param db MetalinksDB database list
#' @return data.frame of sensing relationships
#' @keywords internal
.get_sensing_proteins <- function(db) {
  # Filter for ligand-receptor relationships (type = lr)
  lr_edges <- db$edges[db$edges$type == "lr", ]

  # Merge with protein info
  lr_edges <- merge(lr_edges, db$proteins[, c("uniprot", "gene_symbol", "protein_type")],
    by = "uniprot", all.x = TRUE
  )

  # Merge with metabolite names
  lr_edges <- merge(lr_edges, db$metabolites[, c("hmdb", "metabolite")],
    by = "hmdb", all.x = TRUE
  )

  lr_edges
}

#' Get Transporter Proteins
#'
#' @param db MetalinksDB database list
#' @param direction Character. "out" for secretion, "in" for uptake
#' @return data.frame of transporter relationships
#' @keywords internal
.get_transporters <- function(db, direction = "out") {
  if (!direction %in% c("in", "out")) {
    stop("direction must be 'in' or 'out'")
  }

  # Filter for transporter proteins with specific direction
  trans_edges <- db$edges[!is.na(db$edges$transport_direction) &
    db$edges$transport_direction == direction, ]

  # Merge with protein info
  trans_edges <- merge(trans_edges, db$proteins[, c("uniprot", "gene_symbol", "protein_type")],
    by = "uniprot", all.x = TRUE
  )

  # Merge with metabolite names
  trans_edges <- merge(trans_edges, db$metabolites[, c("hmdb", "metabolite")],
    by = "hmdb", all.x = TRUE
  )

  trans_edges
}

#' Get Extracellular Metabolites
#'
#' @param db MetalinksDB database list
#' @return Character vector of HMDB IDs for extracellular metabolites
#' @keywords internal
.get_extracellular_metabolites <- function(db) {
  unique(db$cell_location$hmdb[db$cell_location$cell_location == "Extracellular"])
}

#' Get Pathway Information for Metabolites
#'
#' @param db MetalinksDB database list
#' @param hmdb_ids Character vector of HMDB IDs
#' @return data.frame of pathway associations
#' @keywords internal
.get_metabolite_pathways <- function(db, hmdb_ids = NULL) {
  if (is.null(hmdb_ids)) {
    return(db$pathway)
  }
  db$pathway[db$pathway$hmdb %in% hmdb_ids, ]
}

#' List Available Metabolites
#'
#' @description Get a list of all metabolites in the database with their properties
#'
#' @param type Character. Filter by interaction type: "all", "signaling" (lr only),
#'        or "metabolic" (pd only)
#'
#' @return data.frame with metabolite information
#' @export
#'
#' @examples
#' \donttest{
#' # Get all metabolites
#' mets <- listMetabolites()
#' head(mets)
#'
#' # Get only signaling metabolites (those with receptors)
#' signaling_mets <- listMetabolites(type = "signaling")
#' }
listMetabolites <- function(type = "all") {
  if (!type %in% c("all", "signaling", "metabolic")) {
    stop("type must be 'all', 'signaling', or 'metabolic'")
  }

  db <- .load_metalinksdb()

  if (type == "all") {
    return(db$metabolites)
  } else if (type == "signaling") {
    lr_hmdb <- unique(db$edges$hmdb[db$edges$type == "lr"])
    return(db$metabolites[db$metabolites$hmdb %in% lr_hmdb, ])
  } else if (type == "metabolic") {
    pd_hmdb <- unique(db$edges$hmdb[db$edges$type == "pd"])
    return(db$metabolites[db$metabolites$hmdb %in% pd_hmdb, ])
  }
}

#' List Available Genes
#'
#' @description Get a list of all genes in the database with their roles
#'
#' @param role Character. Filter by role: "all", "enzyme", "receptor", "transporter"
#'
#' @return data.frame with gene information
#' @export
#'
#' @examples
#' \donttest{
#' # Get all genes
#' genes <- listGenes()
#'
#' # Get only receptor genes
#' receptors <- listGenes(role = "receptor")
#' }
listGenes <- function(role = "all") {
  if (!role %in% c("all", "enzyme", "receptor", "transporter")) {
    stop("role must be 'all', 'enzyme', 'receptor', or 'transporter'")
  }

  db <- .load_metalinksdb()

  proteins <- db$proteins

  # Clean protein_type (remove quotes if present)
  proteins$protein_type_clean <- sapply(proteins$protein_type, .clean_protein_type)

  if (role == "all") {
    return(proteins[, c("uniprot", "gene_symbol", "protein_type")])
  } else if (role == "enzyme") {
    # Enzymes are typically those without a specific receptor/transporter type
    # or explicitly marked as enzyme
    enzyme_match <- proteins$protein_type_clean == "enzyme" | is.na(proteins$protein_type_clean)
    return(proteins[enzyme_match, c("uniprot", "gene_symbol", "protein_type")])
  } else if (role == "receptor") {
    receptor_types <- c("gpcr", "lgic", "nhr", "vgic", "catalytic_receptor", "other_ic")
    receptor_match <- proteins$protein_type_clean %in% receptor_types
    return(proteins[receptor_match, c("uniprot", "gene_symbol", "protein_type")])
  } else if (role == "transporter") {
    transporter_match <- proteins$protein_type_clean == "transporter"
    return(proteins[transporter_match, c("uniprot", "gene_symbol", "protein_type")])
  }
}

#' Get Database Statistics
#' @keywords internal
.get_db_stats <- function(db) {
  list(
    n_metabolites = nrow(db$metabolites),
    n_proteins = nrow(db$proteins),
    n_edges = nrow(db$edges),
    n_lr = sum(db$edges$type == "lr"),
    n_pd = sum(db$edges$type == "pd"),
    n_pathways = length(unique(db$pathway$pathway)),
    n_diseases = length(unique(db$disease$disease))
  )
}
