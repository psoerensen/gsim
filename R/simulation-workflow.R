#' Run the shared standalone simulation workflow from R
#'
#' This additive workflow accepts the same named settings as simulation.yml.
#' It leaves gsim(), gsim_pedigree() and gsim_records() behavior unchanged.
#' Inputs are a parent-first FAM pedigree, one-chromosome BIM, founder phase
#' text (ID and two bit strings), and ordered marker-effect text. The simulation
#' shares native meiosis, packed storage and identity-keyed RNG with standalone.
#' Gaussian model recipes use schema gsim.gaussian_model/v1 with model and
#' design_inputs, or inputs$specification for an already prepared model.
#' Model preparation uses the installed public gsuite command; prepared Gaussian
#' execution uses the same native gsim core without Python or R RNG draws.
#' @param simulation Named list with the selected schema and its settings.
#' @param output New output directory; its parent must exist.
#' @param base_directory Base directory for relative input paths.
#' @return Invisibly the normalized output directory, containing CSV/JSON
#'   and, in marker mode, PLINK files.
#' @export
gsim_simulate <- function(simulation, output, base_directory = ".") {
  if(is.list(simulation)&&identical(simulation$schema,"gsim.gaussian_model/v1"))
    return(.gsim_gaussian_simulate(simulation,output,base_directory))
  allowed <- c("schema","inputs","traits","seed","genotyped","phenotyped","records_per_individual","visits",
               "missing_rate","residual_covariance","maximum_memory_bytes")
  if(!is.list(simulation)||is.null(names(simulation))||anyDuplicated(names(simulation))||
     any(!names(simulation)%in%allowed)||!identical(simulation$schema,"gsim.simulation/v1"))
    stop("Invalid simulation schema or unknown keys.",call.=FALSE)
  if(all(c("records_per_individual","visits") %in% names(simulation)))
    stop("Specify records_per_individual or legacy visits, not both.",call.=FALSE)
  input <- simulation$inputs
  if(!is.list(input)||!setequal(names(input),c("pedigree","variants","founders","effects"))||anyDuplicated(names(input)))
    stop("inputs must name pedigree, variants, founders and effects.",call.=FALSE)
  quote <- function(x) encodeString(enc2utf8(x),quote='"')
  lines <- "GSIM_SIMULATION 1"
  for(k in c("pedigree","variants","founders","effects")) {
    v<-input[[k]];if(!is.character(v)||length(v)!=1L||is.na(v)||!nzchar(v))stop("Input paths must be strings.")
    absolute<-grepl("^(/|[A-Za-z]:[/\\\\])",v)
    path<-normalizePath(if(absolute)v else file.path(base_directory,v),winslash="/",mustWork=TRUE)
    if(dir.exists(path))stop("Input must be a regular file.")
    lines<-c(lines,paste(k,quote(path)))
  }
  traits<-simulation$traits
  if(!is.character(traits)||!length(traits)||anyNA(traits))stop("traits must be strings.")
  if(any(!grepl("^[A-Za-z0-9_.-]+$",traits))||any(traits=="0"))stop("Invalid trait identity.")
  lines<-c(lines,paste("traits",quote(paste(traits,collapse=" "))))
  for(k in c("seed","genotyped","phenotyped","records_per_individual","visits","maximum_memory_bytes"))if(!is.null(simulation[[k]])) {
    v<-simulation[[k]];if(!is.numeric(v)||length(v)!=1||!is.finite(v)||v!=floor(v)||v<0||v>2^53-1)stop(paste(k,"must be an exactly represented nonnegative integer"))
    lines<-c(lines,paste(k,quote(sprintf("%.0f",v))))
  }
  if(!is.null(simulation$missing_rate)) {
    v<-simulation$missing_rate;if(!is.numeric(v)||length(v)!=1||!is.finite(v))stop("missing_rate must be finite numeric")
    lines<-c(lines,paste("missing_rate",quote(sprintf("%.17g",v))))
  }
  if(!is.null(simulation$residual_covariance)) {
    a<-simulation$residual_covariance
    if(!is.matrix(a)||!is.numeric(a)||!identical(dim(a),rep(length(traits),2))||any(!is.finite(a)))stop("Residual covariance must be a trait-ordered numeric matrix.")
    lines<-c(lines,paste("residual_covariance",quote(paste(sprintf("%.17g",as.vector(t(a))),collapse=" "))))
  }
  if(!is.character(output)||length(output)!=1||is.na(output)||!nzchar(output))stop("output must be a path")
  target<-file.path(normalizePath(dirname(output),winslash="/",mustWork=TRUE),basename(output))
  .Call(C_gsim_simulation_workflow,paste0(paste(lines,collapse="\n"),"\n"),target)
  invisible(target)
}
