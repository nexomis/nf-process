# Contribution Rules and Conventions: PROCESS

Objectives:

- Ensure consistency and best practices in process definitions
- Facilitate maintainability and readability of Nextflow processes
- Promote efficient and reproducible pipeline development

## 0. Project Structure

- Each process should accomplish a single task or a closely related set of tasks.

- Store processes in individual files nmaed `main.nf` within a directory structure that reflects the tool and analysis step.

- Use clear, descriptive names for processes
  - Use lower_snake_case for directory holding the process script `main.nf`
  - Optionally use sub-directory for tools with subcommand such as `kraken2` (`index` vs `classify`)
  - following `UPPER_SNAKE_CASE` convention for process name in `main.nf` (e.g., `FASTP`, `KRAKEN2`).

Example structure:
```
process/
  fastp/
    main.nf
  kraken2/
    main.nf
    build/
      main.nf
  spades/
    main.nf
```

Example `main.nf`:
```
process FASTP {
  ...
}
```

## 1. Input and Output Handling

### 1.1 Input declaration

- Use `tuple val(meta), path(...)` only for sample-specific inputs.
- Specify input `arity` and use `stageAs` to avoid conflicts.

Example: 
```
  input:
  tuple val(meta), path(reads, arity: 1..2, stageAs: 'input_raw/*')
```

- When comparing to a reference or database, separate inputs shall not use the same `meta`.
Example:
```
  input:
  tuple val(meta), path(reads, arity: 1..2, stageAs: 'input_raw/*')
  tuple val(meta2), path(index, arity: 1, stageAs: 'input_ref/*')
```

### 1.2 output declaration

- Use `tuple val(meta), path(...)` for sample-specific outputs.
- Include, if possible, `meta.id` in file names or file path to facilitate sample to file association.
- Do not rename or move files to change its name as it's not compatible with S3-based storage.
- When possible and for big files, choose a compressed format output. 

Example: 
```
output:
tuple val(meta), path("${meta.id}*.fq.gz", arity: 1..2), emit: reads
tuple val(meta), path("${meta.id}.json", arity: 1), emit: report
```

- Files with index are defined as flatten with meta: 

Example:
```
  output:
  tuple val(meta), path("${meta.id}.bam"), path("${meta.id}.bam.bai"), emit: bam
```

## 2. Metadata map `meta`

- Use meta map for sample-specific information.
- Include at least meta.id as a unique identifier.

### 2.1 Attributes names by convention

- `meta.read_type`:
  - `SR` for single-read
  - `PE` for pair-end 
  - `LR` for long-read 
  - `spring` for compressed SR or PE  

- `meta.strand`: Type of stranded library for RNA-Seq

- `meta.ref_id`: Reference identifier.

- `meta.args_{{tool}}`: Tool/process-specific arguments.
  - `meta.strand` - in kallisto format: fr-stranded, rf-stranded, or unstranded.
  - `meta.is_3prime` - boolean relative to library prep (false assumes 'full_length' sequencing. Note: '5prime' libraries are not yet supported).
  - `meta.frag_size` and `meta.frag_size_sd` - float relative to mean/sd fragment size.


### 2.2 To avoid with `meta`

- Do not embed file in meta.
- Only meta.id should be mandatory in process definitions; other attributes should have fallbacks.

## 3 Script & command line

### 3.1 Pipeline-specific arguments

- Define default arguments within the process.
- Use `task.ext.args` for user-overridable parameters.
- Provide fallback to empty string or defaults

Example:
```
  script:
  def default_args = "--trim_poly_g --cut_right"
  """
  fastp --thread $task.cpus ${task.ext.args ?: default_args} \\
    --json ${meta.id}.json \\
    $in_args $out_args
```

### 3.2 Sample-Specific Arguments

- Use `meta.args_{{tool}}` for sample-specific arguments.
- Provide fallbacks for when these are not defined.

Example:
```
  spades.py --threads ${task.cpus} \\
    -o ${meta.id} \\
    ${meta.args_spades ?: ''} \\
    ${task.ext.args ?: ''} \\
```

### 3.3 Command

- Avoid capturing `stderr` and `stdout` unless necessary for output generation.
- Keep complex Groovy logic separate from the bash script.

## 4. Ressource management

- Specify CPU and memory requirements when needed using labels such as cpu_med, mem_8G.
- see https://github.com/nexomis/nf-config/blob/main/process/labels.config
- Override default requirements through pipeline-specific configuration if necessary.

Example:
```
  label 'cpu_med'
  label 'mem_8G'
  script:
  def memory_in_mb = task.memory.toMega()
  """
  fastqc --threads \$task.cpus --memory \$memory_in_mb \$reads
  """
```

## 5. Containerization

- Specify container images with version-specific tags for reproducibility.
- Prefer lightweight containers, and consider creating new images if necessary.
- Use `params.biocontainers_registry` to define the biocontainer registry (aws or quay.io) with a fallback

Example:
```
  container "${params.biocontainers_registry ?: 'quay.io'}/biocontainers/kallisto:0.50.1--h6de1650_2"
```

## 6. Misc

- Provide as many relevant outputs as possible.
- Avoid Internal Logic, keep process execution logic in the workflow, not within the process.
- Dot not define publishDir directive at the process level. It needs to be defined at the pipeline level.
- Include stubs for workflow developments.

### 6.1 Documentation

- Provide brief descriptions of process purpose, inputs, and outputs in comments.
- Document any non-obvious behavior or requirements

### 6.2 Global rules and patterns

- Take in account the global rules defined here: https://github.com/nexomis/nf-template/blob/main/CONTRIBUTE.md#typography
- Use empty files to manage optional inputs.

Example: 
```
process QUAST {
  ...
  input:
  tuple val(meta), path(assembly, stageAs: "inputs/assembly.fa")
  tuple val(meta2), path(ref_fa, stageAs: "inputs/reference.fa")
  tuple val(meta3), path(bam, stageAs: "inputs/aln.bam"), path(bai, stageAs: "inputs/aln.bam.bai") 
  ...
  script:
  def args_ref = ref_fa.size() > 1 ? "-r inputs/reference.fa" : ""
  def args_bam = bam.size() > 1 ? "--bam inputs/aln.bam" : ""
  """
  quast.py \\
    --output-dir ${meta.id} \\
    --labels ${meta.id} \\
    --threads $task.cpus \\
    $args_ref \\
    $args_bam \\
    ${task.ext.args ?: ''} \\
    $assembly \\
    2> ${meta.id}.log
  """
}

```

