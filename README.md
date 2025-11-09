# Contribution Rules and Conventions: PROCESS

Objectives:

- Ensure consistency and best practices in process definitions
- Facilitate maintainability and readability of Nextflow processes
- Promote efficient and reproducible pipeline development

## 0. Project Structure

- Each process should accomplish a single task or a closely related set of tasks.

- Store processes in individual files named `main.nf` within a directory structure that reflects the tool and analysis step.

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
- Prefer `meta.label ?: meta.id` in file names for user-friendly output naming.
- Use `meta.id` for channel operations where uniqueness is critical.
- **Do not** rename or move files to change their names, as this is poorly supported with S3-based storage.
- When possible, and especially for large files, output in compressed formats.

**Example:** 
```
output:
tuple val(meta), path("${meta.label ?: meta.id}*.fq.gz", arity: 1..2), emit: reads
tuple val(meta), path("${meta.label ?: meta.id}.json", arity: 1), emit: report
```

**Note:** While `meta.label ?: meta.id` is preferred for file names, continue using `meta.id` for channel operations where uniqueness is critical.

- For files with indexes, define them as flattened with `meta`. 

**Example:**
```
  output:
  tuple val(meta), path("${meta.label ?: meta.id}.bam"), path("${meta.label ?: meta.id}.bam.bai"), emit: bam
```

## 2. Metadata Map `meta`

The `meta` map carries sample-specific information through the pipeline. It serves as the primary mechanism for tracking sample identity, configuration, and metadata throughout the workflow.

> **Note:** See section 2.2 for a complete list of reserved meta attributes to avoid naming conflicts.

### 2.1 Core Principles & Usage

#### `meta.id` - Unique Identifier

- **MUST be unique** across all samples in the pipeline
- Used for **channel operations** such as joins, grouping, and combinations
- Critical for dataflow operations where sample matching is required
- Always use `meta.id` in channel operations to ensure correct sample tracking

**Example:**
```groovy
// Channel join relies on matching meta.id
ch_reads
  .map{ it -> [it[0].id, it]}
  .join(
    ch_qc_reports
      .map{ it -> [it[0].id, it]},
    by: 0
   )
   .map{ it -> [it[1], it[2]] } // Joins on first element (meta.id)
   // or better: .map{ it -> [it[1][0], it[1][1], it[2][1]] }
```

#### `meta.label` - Display Name

- **Preferred name** for output files and reports
- Provides a user-friendly alternative to `meta.id` for file naming
- **Should be unique if provided** to avoid file name conflicts in published outputs
- Use the pattern `${meta.label ?: meta.id}` for output file names (ensures uniqueness via fallback to `meta.id`)
- Use the pattern `tag "${meta.label ?: meta.id}"` in process definitions

**Example:**
```groovy
process ANALYSIS {
    tag "${meta.label ?: meta.id}"
    
    output:
    tuple val(meta), path("${meta.label ?: meta.id}.report.txt"), emit: report
}
```

**Important:** 
- The pattern `${meta.label ?: meta.id}` ensures uniqueness because it falls back to the unique `meta.id` when `meta.label` is not provided
- If users provide `meta.label`, they should ensure it is also unique across samples to prevent file overwriting
- Continue using `meta.id` for channel operations (joins, grouping) where uniqueness is critical

#### Mandatory vs Optional Attributes

- Only `meta.id` should be mandatory in process definitions
- **All other attributes must have fallbacks** using the elvis operator (`?:`)
- This ensures processes work correctly even when optional metadata is not provided

### 2.2 Reserved Meta Attributes

The following meta attributes have standardized meanings and should be used consistently across processes. Avoid using these names for other purposes.

#### Core Identifiers

- **`meta.id`** - Unique identifier (mandatory)
- **`meta.label`** - Display name for outputs and reports (optional, falls back to `meta.id`)

#### Read/Data Type Information

- **`meta.read_type`** - Read type indicator:
  - `SR` for single-read
  - `PE` for pair-end
  - `LR` for long-read
  - `unknown` (e.g. for compressed/sra)
- **`meta.strand`** - Strandedness (kallisto format):
  - `fr-stranded` - forward-reverse stranded
  - `rf-stranded` - reverse-forward stranded
  - `unstranded` - unstranded
- **`meta.is_3prime`** - Boolean indicating 3' library prep
- **`meta.frag_size`** - Mean fragment size (float)
- **`meta.frag_size_sd`** - Fragment size standard deviation (float)

#### Reference & Database

- **`meta.ref_id`** - Reference genome/database identifier

#### Tool-Specific Arguments

> **These attributes are not mandatory** however will be useful. note that they can co-exist with `task.ext.args`.

- **`meta.args_<processname>`** - Process-specific command-line arguments
  - Example: `meta.args_fastp` for the FASTP process
  - Example: `meta.args_kraken2` for the KRAKEN2 process
- **`meta.args_<instance_name>`** - Arguments for aliased process instances
  - Example: `meta.args_fastp_trim` for `FASTP as FASTP_TRIM`
  - Example: `meta.args_fastp_dedup` for `FASTP as FASTP_DEDUP`

> **NOTE** Attention, `<instance_name>` shall not be used in process script. Workflow developers shall use `map` to pass it to the `<process_name>` key before process input.

#### Optional Resource Hints (NOT DEFAULT)

⚠️ **These attributes are NOT the recommended default approach.** Use only when users need fine-grained per-sample resource control.

- **`meta._cpus`** - Per-sample CPU override (integer)
- **`meta._memory`** - Per-sample memory override (MemoryUnit, e.g., `32.GB`)
- **`meta._time`** - Per-sample time limit override (Duration, e.g., `4.hour`)

**Note:** See Section 4.2, Pattern 6 for details on when and how to use resource hints.

#### Other exemples (not reserved) Sample-Specific Attributes

- **`meta.kmer_size`** - K-mer size for assembly processes
- **`meta.genome_size`** - Genome size category (e.g., `large`, `medium`, `small`)

### 2.3 Restrictions with `meta`

- **Do not** embed files in `meta` (use separate path inputs instead)
- Only `meta.id` should be mandatory in process definitions; other attributes must have fallbacks
- Always clone meta in `map` operations when modifying it.
- **Avoid naming conflicts** with reserved attributes listed in section 2.2
- When introducing new meta attributes, document them and consider if they should be added to the reserved list

## 3 Scripts & Command Line Tools

### 3.1 Process-Specific Arguments (NF 25.10+)

**Important:** In Nextflow 25.10+ with strict syntax, process arguments come from **two sources**:

1. **params-dependent arguments** (via `meta.args_<processname>`):
   - Set in entry workflow by mapping params to meta
   - Must be passed explicitly to processes
   - Allows different arguments for different workflow contexts
   
2. **params-independent arguments** (via `task.ext.*`):
   - Set in config files (`conf/ext.config`)
   - Apply regardless of workflow parameters
   - Used for static/default behavior

**Pattern for meta.args_\<processname\>:**

The standardized naming convention is:
- Process `FASTP` reads from `meta.args_fastp`
- Process `KRAKEN2` reads from `meta.args_kraken2`
- Process `SPADES` reads from `meta.args_spades`

**Example in entry workflow (main.nf):**
```groovy
// Map workflow params to meta
ch_samples
    .map { meta, reads ->
        def new_meta = meta.clone()
        new_meta.args_fastp = params.fastp_opts  // Will be read by FASTP process
        new_meta.args_kraken2 = params.kraken2_opts
        [new_meta, reads]
    }
    .set { ch_input }
```

**Example in process:**
```groovy
process FASTP {
    tag "${meta.label ?: meta.id}"
    
    input:
    tuple val(meta), path(reads)
    
    output:
    tuple val(meta), path("${meta.label ?: meta.id}.fq.gz"), emit: reads
    tuple val(meta), path("${meta.label ?: meta.id}.json"), emit: json
    
    script:
    // Workflow-dependent args from meta
    def workflow_args = meta.args_fastp ?: ''
    // Workflow-independent defaults from config
    def default_args = task.ext.default_args ?: '--trim_poly_g --cut_right'
    
    """
    fastp --thread ${task.cpus} \\
        ${default_args} \\
        ${workflow_args} \\
        --json ${meta.label ?: meta.id}.json \\
        -i ${reads[0]} -o ${meta.label ?: meta.id}.fq.gz
    """
}
```

### 3.2 Handling Multiple Process Instances

When a process is called multiple times in a workflow with different names, use meta keys matching the process instance names:

**Process instances in workflow/subworkflow:**
```groovy
include { FASTP as FASTP_TRIM } from './process/fastp/main.nf'
include { FASTP as FASTP_DEDUP } from './process/fastp/main.nf'
```

**Entry workflow:**
```groovy
ch_samples
    .map { meta, reads ->
        def new_meta = meta.clone()
        new_meta.args_fastp_trim = params.trim_opts      // For FASTP_TRIM
        new_meta.args_fastp_dedup = params.dedup_opts    // For FASTP_DEDUP
        [new_meta, reads]
    }
    .set { ch_input }
```

**Subworkflow:**
```groovy
// First call: map args_fastp_trim → args_fastp
samples
    .map { meta, reads ->
        def new_meta = meta.clone()
        new_meta.args_fastp = meta.args_fastp_trim ?: ''
        [new_meta, reads]
    }
    .set { ch_for_trim }

FASTP_TRIM(ch_for_trim)

// Second call: map args_fastp_dedup → args_fastp
FASTP_TRIM.out.reads
    .map { meta, reads ->
        def new_meta = meta.clone()
        new_meta.args_fastp = meta.args_fastp_dedup ?: ''
        [new_meta, reads]
    }
    .set { ch_for_dedup }

FASTP_DEDUP(ch_for_dedup)
```

**Key principle:** Meta keys follow the pattern `meta.args_<instance_name_lowercase>` where instance name is the process alias used in the workflow.

### 3.3 Sample-Specific Arguments (Other Meta Attributes)

Use other meta attributes for sample-specific behavior:

**Example:**
```groovy
process SPADES {
    tag "${meta.label ?: meta.id}"
    
    input:
    tuple val(meta), path(reads)
    
    output:
    tuple val(meta), path("${meta.label ?: meta.id}/*"), emit: assembly
    
    script:
    def args = meta.args_spades ?: ''
    def kmer_size = meta.kmer_size ? "-k ${meta.kmer_size}" : ''
    
    """
    spades.py --threads ${task.cpus} \\
        -o ${meta.label ?: meta.id} \\
        ${args} \\
        ${kmer_size} \\
        ${task.ext.args ?: ''}
    """
}
```

### 3.3 Script Best Practices

- **Avoid** capturing `stderr` and `stdout` unless necessary for output generation.
- Keep complex Groovy logic separate from the Bash script to maintain readability.

## 4. Resource Management

### 4.1 Static Resource Allocation

- Override default requirements through pipeline-specific configuration if necessary [conf/resources.config](https://github.com/nexomis/nf-template/blob/main/conf/resources.config).

**Example:**
```groovy
process FASTQC {
  cpus 8
  memory 8.GB
  
  script:
  def memory_in_mb = task.memory.toMega()
  """
  fastqc --threads ${task.cpus} --memory ${memory_in_mb} ${reads}
  """
}
```

### 4.2 Dynamic Resource Allocation

> **⚠️ Important: Not the Default Approach**
>
> The strategies described in this section are **NOT the recommended default** for resource allocation. Static resource allocation via configuration files (`conf/resources.config`) is the baseline approach for most pipelines.
>
> Dynamic allocation should be used **selectively** for specific scenarios where static allocation is insufficient. Each pattern below indicates when it is appropriate to use.

Dynamic directives allow process resources to be adjusted based on task attempt, input characteristics, or metadata. This is useful when different instances of the same process have very different resource requirements.

**Important (NF 25.10+):** With strict syntax, dynamic directives in processes no longer require closures. Use direct expressions:

```groovy
// v25.10.0 strict syntax
memory 32.GB * task.attempt
time 4.hour * task.attempt
cpus meta.read_type == 'LR' ? 16 : 8
```

#### Pattern 1: Retry with Increasing Resources

Automatically increase resources when a task fails and is re-executed:

```groovy
process ASSEMBLY {
    memory 32.GB * task.attempt
    time 4.hour * task.attempt
    
    maxRetries 3
    
    input:
    tuple val(meta), path(reads)
    
    script:
    """
    assembly_tool --threads ${task.cpus} --memory ${task.memory.toGiga()} ${reads}
    """
}
```

**How it works:**
- First attempt: `task.attempt = 1` → 32 GB memory, 4 hours
- Second attempt: `task.attempt = 2` → 64 GB memory, 8 hours  
- Third attempt: `task.attempt = 3` → 96 GB memory, 12 hours

#### Pattern 2: Resources Based on Input File Size

Adjust resources based on the size of input files:

```groovy
process ALIGNMENT {
    memory 8.GB + 1.GB * Math.ceil(reads.size() / 1024 ** 3)
    
    input:
    tuple val(meta), path(reads)
    
    script:
    """
    aligner --threads ${task.cpus} ${reads}
    """
}
```

**How it works:**
- Base memory: 8 GB
- Additional memory: 1 GB per GB of input file size (rounded up)
- Example: 5.2 GB input file → 8 + 6 = 14 GB memory requested

#### Pattern 3: Resources Based on Meta Attributes

Use metadata to determine resource requirements:

```groovy
process MAPPING {
    cpus meta.read_type == 'LR' ? 16 : 8
    memory meta.genome_size == 'large' ? 64.GB : 32.GB
    
    input:
    tuple val(meta), path(reads), path(index)
    
    script:
    """
    mapper --threads ${task.cpus} ${reads} ${index}
    """
}
```

#### Pattern 4: Dynamic Resources with Previous Execution Trace (NF 24.10+)

Adjust resources based on metrics from previous task attempts:

```groovy
process VARIANT_CALLING {
    memory task.attempt > 1 ? task.previousTrace.memory * 1.5 : 16.GB
    time task.attempt > 1 ? task.previousTrace.realtime * 2 : 2.hour
    
    maxRetries 3
    
    input:
    tuple val(meta), path(bam), path(bai)
    
    script:
    """
    variant_caller --threads ${task.cpus} ${bam}
    """
}
```

**How it works:**
- First attempt: 16 GB memory, 2 hours
- Subsequent attempts: 1.5× previous memory, 2× previous runtime
- More efficient than fixed multipliers as it adapts to actual resource usage

#### Pattern 5: Combined Dynamic Resources

Combine multiple strategies for robust resource allocation:

```groovy
process DENOVO_ASSEMBLY {
    // Base memory on input size, increase on retry
    memory (16.GB + 2.GB * Math.ceil(reads.size() / 1024 ** 3)) * task.attempt
    
    // CPUs based on read type
    cpus meta.read_type == 'LR' ? 32 : 16
    
    // Time based on genome size category
    time (meta.genome_size == 'large' ? 24.hour : 8.hour) * task.attempt
    
    maxRetries 2
    
    input:
    tuple val(meta), path(reads)
    
    script:
    """
    assembler --threads ${task.cpus} --memory ${task.memory.toGiga()} ${reads}
    """
}
```

**When to use:**
- Complex workflows with multiple varying factors
- Processes where both input characteristics and retry strategies are needed

#### Pattern 6: Per-Sample Resource Hints via Meta (Advanced)

> **⚠️ NOT Recommended as Default**
>
> This pattern should **only** be used when users need fine-grained, per-sample resource control (e.g., known problematic samples requiring specific resources). This is an **advanced use case** and should not be the standard approach.

Allow users to specify resources per sample through the samplesheet/metadata using reserved meta attributes:

```groovy
process CUSTOM_ANALYSIS {
    tag "${meta.label ?: meta.id}"
    
    // Optional: Allow meta to override defaults
    cpus meta._cpus ?: 8
    memory meta._memory ?: 16.GB
    time meta._time ?: 2.hour
    
    input:
    tuple val(meta), path(input_file)
    
    output:
    tuple val(meta), path("${meta.label ?: meta.id}.result.txt"), emit: results
    
    script:
    """
    analysis_tool --threads ${task.cpus} \\
        --memory ${task.memory.toGiga()} \\
        ${input_file} > ${meta.label ?: meta.id}.result.txt
    """
}
```

**When to use this pattern:**
- Users have prior knowledge of specific samples needing more resources
- Resource requirements cannot be predicted from file size or metadata categories
- Fine-grained control is explicitly required by the workflow design
- **NOT for general-purpose pipelines** - prefer static config or other dynamic patterns

**When NOT to use:**
- As the default resource allocation strategy
- When file size or metadata categories (Pattern 2 or 3) would suffice
- When retry strategies (Pattern 1 or 4) can handle resource failures
- For public/shared pipelines where users should not need to specify resources per sample

**Naming convention (see Section 2.2):**
- `meta._cpus` - CPU count (integer)
- `meta._memory` - Memory amount (MemoryUnit, e.g., `32.GB`)
- `meta._time` - Time limit (Duration, e.g., `4.hour`)

**Usage in entry workflow:**
```groovy
// In entry workflow when creating meta from samplesheet or params
ch_samples
    .map { meta, reads ->
        def new_meta = meta.clone()
        // Optional resource hints for specific samples
        // Can come from samplesheet columns or params
        if (params.sample_resources?.containsKey(meta.id)) {
            new_meta._cpus = params.sample_resources[meta.id].cpus
            new_meta._memory = params.sample_resources[meta.id].memory
        }
        [new_meta, reads]
    }
```

**Important:** Always provide sensible defaults using the elvis operator (`?:`) to ensure the process works correctly even when these attributes are not set.

#### Best Practices for Dynamic Resources

**✅ DO:**
- Use `task.attempt` for retry strategies
- Calculate resources based on input file sizes for predictable scaling
- Use `meta` attributes to categorize different resource profiles
- Set reasonable `maxRetries` to avoid infinite retry loops
- Consider `task.previousTrace` for adaptive resource scaling (NF 24.10+)

**❌ DON'T:**
- Don't use params directly in process directives (breaks modularity)
- Don't set unlimited retries without maxRetries
- Don't use very aggressive multipliers (e.g., `10 * task.attempt`)

#### Exit Codes Reference

Common exit codes for resource-related failures:
- `137` - SIGKILL (often out of memory - OOM)
- `138` - Container killed (resource limit exceeded)
- `139` - SIGSEGV (segmentation fault, may indicate memory issues)
- `140` - Child process killed

When implementing custom error strategies, these codes can be used to detect resource-related failures:

```groovy
errorStrategy task.exitStatus in 137..140 ? 'retry' : 'terminate'
```

#### Example: Complete Process with Dynamic Resources

```groovy
process KRAKEN2_CLASSIFY {
    container 'ghcr.io/nexomis/kraken2:2.1.3'
    tag "${meta.label ?: meta.id}"
    
    cpus 16
    
    // Dynamic memory based on database size and retry attempt
    memory 1.GB * Math.ceil(kraken2_db.size() / 1024 ** 3) * 1.2 * task.attempt
    
    time 1.hour * task.attempt
    
    maxRetries 2
    
    input:
    tuple val(meta), path(reads, arity: 1..2)
    tuple val(meta2), path(kraken2_db)
    
    output:
    tuple val(meta), path("${meta.label ?: meta.id}.kraken2.txt"), emit: report
    tuple val(meta), path("${meta.label ?: meta.id}_unclassified*.fq.gz", arity: 1..2), emit: unclassified
    
    script:
    def is_paired = reads.size() == 2
    def input_arg = is_paired ? "--paired ${reads[0]} ${reads[1]}" : "${reads[0]}"
    def unclass_out = is_paired ? 
        "--unclassified-out ${meta.label ?: meta.id}_unclassified#.fq" :
        "--unclassified-out ${meta.label ?: meta.id}_unclassified.fq"
    
    """
    kraken2 \\
        --db ${kraken2_db} \\
        --threads ${task.cpus} \\
        --output ${meta.label ?: meta.id}.kraken2.txt \\
        ${unclass_out} \\
        ${input_arg}
    
    gzip ${meta.label ?: meta.id}_unclassified*.fq
    """
}
```

## 5. Containerization

- Specify container images with **version-specific tags** for reproducibility.
- Prefer lightweight containers; consider creating new images if necessary.

**Example:**
```
  container "quay.io/biocontainers/kallisto:0.50.1--h6de1650_2"
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
  tag "${meta.label ?: meta.id}"
  
  input:
  tuple val(meta), path(assembly, stageAs: "inputs/assembly.fa")
  tuple val(meta2), path(ref_fa, stageAs: "inputs/reference.fa")
  tuple val(meta3), path(bam, stageAs: "inputs/aln.bam"), path(bai, stageAs: "inputs/aln.bam.bai") 
  
  output:
  tuple val(meta), path("${meta.label ?: meta.id}/*"), emit: results
  tuple val(meta), path("${meta.label ?: meta.id}.log"), emit: log
  
  script:
  def args_ref = ref_fa.size() > 1 ? "-r inputs/reference.fa" : ""
  def args_bam = bam.size() > 1 ? "--bam inputs/aln.bam" : ""
  """
  quast.py \\
    --output-dir ${meta.label ?: meta.id} \\
    --labels ${meta.label ?: meta.id} \\
    --threads $task.cpus \\
    $args_ref \\
    $args_bam \\
    ${task.ext.args ?: ''} \\
    $assembly \\
    2> ${meta.label ?: meta.id}.log
  """
}

```
