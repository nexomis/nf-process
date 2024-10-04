# Contribution Rules and Conventions: PROCESS

Objectives:

- Ensure consistency and best practices in process definitions
- Facilitate maintainability and readability of Nextflow processes
- Promote efficient and reproducible pipeline development

## 0. Project Structure

- Each process should accomplish a single task or a closely related set of tasks.

- Store processes in individual files nmaed `main.nf` within a directory structure that reflects the tool and analysis step.

- Use clear, descriptive names for processes:
  - Use **lower_snake_case** for directory holding the process script `main.nf`
  - Optionally use sub-directory for tools with subcommand such as `kraken2` (`index` vs `classify`)
  - following **UPPER_SNAKE_CASE** convention for process name in `main.nf` (e.g., `FASTP`, `KRAKEN2`).

**Example structure:**
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

**Example** `main.nf`:
```
process FASTP {
  ...
}
```

## 1. Input and Output Handling

### 1.1 Input declaration

- Use `tuple val(meta), path(...)` only for sample-specific inputs.
- Specify input `arity` and use `stageAs` to avoid conflicts.

**Example:** 
```
  input:
  tuple val(meta), path(reads, arity: 1..2, stageAs: 'input_raw/*')
```

- When comparing to a reference or database, separate inputs shall not use the same `meta`.

**Example:**
```
  input:
  tuple val(meta), path(reads, arity: 1..2, stageAs: 'input_raw/*')
  tuple val(meta2), path(index, arity: 1, stageAs: 'input_ref/*')
```

### 1.2 Output declaration

- Use `tuple val(meta), path(...)` for sample-specific outputs.
- Include `meta.id` in file names or file path to facilitate sample-to-file association.
- **Do not** rename or move files to change their names, as this is poorly supported with S3-based storage.
- When possible, and especially for large files, output in compressed formats.

**Example:** 
```
output:
tuple val(meta), path("${meta.id}*.fq.gz", arity: 1..2), emit: reads
tuple val(meta), path("${meta.id}.json", arity: 1), emit: report
```

- For files with indexes, define them as flattened with `meta`. 

**Example:**
```
  output:
  tuple val(meta), path("${meta.id}.bam"), path("${meta.id}.bam.bai"), emit: bam
```

## 2. Metadata Map `meta`

- Use the `meta` map for sample-specific information.
- Include at least `meta.id` as a unique identifier.

### 2.1 Some meta attributes

Below are some example of names that can be used.

- `meta.read_type`:
  - `SR` for single-read
  - `PE` for pair-end 
  - `LR` for long-read 
  - `spring` for compressed SR or PE  
- `meta.ref_id`: Reference identifier.
- `meta.args_{{tool}}`: Tool/process-specific arguments.
- `meta.strand` - in kallisto format: fr-stranded, rf-stranded, or unstranded.
- `meta.is_3prime` - boolean relative to library prep.
- `meta.frag_size` and `meta.frag_size_sd` - float relative to mean/sd fragment size.


### 2.2 Restrictions with `meta`

- **Do not** embed file in `meta`.
- Only `meta.id` should be mandatory in process definitions; other attributes must have fallbacks.

## 3 Scripts & Command Line Tools

### 3.1 Pipeline-specific arguments

- Define default arguments within the process.
- Use `task.ext.args` for user-overridable parameters.
- Provide fallback to empty string or defaults

**Example:**
```
  script:
  def default_args = "--trim_poly_g --cut_right"
  """
  fastp --thread $task.cpus ${task.ext.args ?: default_args} \\
    --json ${meta.id}.json \\
    $in_args $out_args
```

### 3.2 Sample-Specific Arguments

- Use `meta.args_<tool>` for sample-specific arguments.
- Provide fallbacks for when these are not defined.

**Example:**
```
  spades.py --threads ${task.cpus} \\
    -o ${meta.id} \\
    ${meta.args_spades ?: ''} \\
    ${task.ext.args ?: ''} \\
```

### 3.3 Script Best Practices

- **Avoid** capturing `stderr` and `stdout` unless necessary for output generation.
- Keep complex Groovy logic separate from the Bash script to maintain readability.

## 4. Ressource management

- Specify CPU and memory requirements using labels such as `cpu_med`, `mem_8G`.
- Refer to [process labels configuration](https://github.com/nexomis/nf-config/blob/main/process/labels.config) fo predefined labels.
- Override default requirements through pipeline-specific configuration if necessary [conf/resources.config](https://github.com/nexomis/nf-template/blob/main/conf/resources.config).

**Example:**
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

- Specify container images with **version-specific tags** for reproducibility.
- Prefer lightweight containers; consider creating new images if necessary.
- Use `params.biocontainers_registry` to define the biocontainer registry (e.g., AWS or Quay.io) with a fallback option.

**Example:**
```
  container "${params.biocontainers_registry ?: 'quay.io'}/biocontainers/kallisto:0.50.1--h6de1650_2"
```

## 6. Miscellaneous Guidelines

- Provide as many relevant outputs as possible to facilitate downstream analysis.
- **Avoid** internal logic within processes; keep execution logic in the workflow definitions.
- **Do not** define the `publishDir` directive at the process level; it should be defined at the pipeline level.
- Include stubs for workflow development.

### 6.1 Documentation

- Provide brief descriptions of the process's purpose, inputs, and outputs in comments.
- Document any non-obvious behavior or special requirements.

### 6.2 Global rules and patterns

- Adhere to the global rules defined in the [typography conventions](https://github.com/nexomis/nf-template/blob/main/CONTRIBUTE.md#typography).
- Use empty files to manage optional inputs.

**Example:**
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
