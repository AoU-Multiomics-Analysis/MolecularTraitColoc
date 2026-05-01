version 1.0

workflow PairwiseColocByGene {
  input {
    String gene_id
    File? eqtl_susie
    File? sqtl_susie
    File? pqtl_susie
  }

  call RunColoc {
    input:
      gene_id = gene_id,
      eqtl_susie = eqtl_susie,
      sqtl_susie = sqtl_susie,
      pqtl_susie = pqtl_susie
  }

  output {
    File coloc_results = RunColoc.results
  }
}


task RunColoc {
  input {
    String gene_id

    File? eqtl_susie
    File? sqtl_susie
    File? pqtl_susie
    String output_prefix 
  }

  command <<<
    Rscript run_pairwise_coloc.R \
      --gene_id "~{gene_id}" \
      --output "~{output_prefix}.pairwise_coloc.tsv" \
      ~{if defined(eqtl_susie) then "--eqtl_susie " + eqtl_susie else ""} \
      ~{if defined(sqtl_susie) then "--sqtl_susie " + sqtl_susie else ""} \
      ~{if defined(pqtl_susie) then "--pqtl_susie " + pqtl_susie else ""}
  >>>

  output {
    File results = "~{output_prefix}.pairwise_coloc.tsv"
  }

  runtime {
    docker: "" 
    cpu: 1
    memory: "8G"
    disks: "local-disk 50 HDD"
  }
}
