# Contribution Rules and Conventions: PROCESS

## 0. Structure
To ensure maximum portability, processes should preferably accomplish single tasks (within reason, for example: `samtools`: `converts`, `sort`, and `index`, or `sam2bam` in an alignment process with `bowtie2`).
Each process should be in an individual file located in a directory named after the tool, a subdirectory named after the analysis step (if relevant), and, if necessary, a subdirectory for specific adaptations (e.g., different versions or Docker containers, different argument considerations, etc.).

## 1. Handling Input and Output File Paths

### 1.1 Input path
Specify the existing state and arity of the input paths. Use the `saveAs` operator to avoid conflicts with process output files. 

### 1.2 No paths in `meta`
Metadata (`meta`) should not contain file paths; it's better to declare them as `path` from the outset.

## 2. File Handling

### 2.1 `mv`/`rename`
Never use `mv` or equivalent: it can produce strange behaviors in S3, especially with large files.

### 2.2 `cp`
Minimize copying (`cp`), especially for large files. To place all outputs of a tool in a specific directory when the tool necessarily generates output in the current directory, you can create the folder and move into it before running the tool.

### 2.3 If usage of `mv`/`cp` is recquired
Should there be cases where the use of `mv`/`cp` is unavoidable on large files, in order to make these worflows usable on the cloud, we will activate `scratch` for the process concerned.

## 3. Handling Paths as Queue Channels

*Unless in exceptional cases,* paths in queue channels should be passed and retrieved as tuples: `val(meta), path(file/dir)`.

### 3.2 Structure of Queue Channels for Paths (input **and** output)
- The first element of these queue channels should be a map named `meta`, which includes at least a unique identifier (`meta.id`).
- The second element is a list of associated file paths.
- Ideally, the format is limited to what has been described (`val(meta), path(file/dir)`), and if multiple paths are needed for the same task (e.g., reads and a reference genome), they can be provided to the process as two distinct input queue channels that follow the described format and are ordered in the same sequence based on their common `meta.id`.
**In some cases, it might be necessary to add additional elements, which is allowed**, (e.g `val(meta), [path(reads_R1, reads_R2)], [kraken_db], [genome_fasta, genome_index]`)

### 3.3 Specific Attributes

#### 3.3.1 Reads
- `meta.read_type` should be defined, with possible values: `spring`, `SR`, or `PE`.
- Other attributes may be useful, and if needed, they should follow these conventions:
  - `meta.strand` - in kallisto format: fr-stranded, rf-stranded, or unstranded.
  - `meta.is_3prime` - boolean relative to library prep (false assumes 'full_length' sequencing. Note: '5prime' libraries are not yet supported).
  - `meta.frag_size` and `meta.frag_size_sd` - float relative to mean/sd fragment size.

Note: When appropriate and not detrimental to the process, it is preferable to provide one file per pair of reads (rather than files per lane or interleaved, for example) with file name suffixes '_1'/'_2' or '_R1'/'_R2'.

## 4. Arguments and Argument Management

### 4.1 Default Arguments
The default values of arguments used in the processes should be defined at the process level and not upstream in the workflow (e.g., `fastp` process).

### 4.2 Sample-Specific Arguments
For sample-specific arguments defined by elements external to the process (e.g., defined by the user in the sample sheet), the arguments will be stored in `meta` according to two possible options: `meta.args` or `meta.args_<toolsName>`.
In any case, all processes should still function if `meta.args` or `meta.args_<toolsName>` are not defined, by providing default values in the process if necessary.
Note: For sample-specific arguments that can be defined within the process, continue to handle this case within the process (e.g., SPAdes input reads: `(reads.size() == 1) ? "-s ${reads}" : "-1 ${reads[0]} -2 ${reads[1]}"`).

#### 4.2.1 `meta.args`
This needs to be updated before calling each subworkflow or process, otherwise, it could lead to execution errors or, worse, incorrect executions (ensure to set it to `null` if nothing is to be added because by default, processes may include `${meta.args ?: ''}` in their command lines).

#### 4.2.2 `meta.args_<toolsName>`
Be cautious if the same tool (e.g., `samtools`, `kraken2`, ...) is used multiple times in the workflow.

## 5. Structure and Portability of Processes

### 5.1 Stub Category
All processes must have a `stub` category that at least creates the necessary elements for a standard workflow.

### 5.2 Process Portability
The logic of whether or not a process should be executed must be defined at the workflow level, not within the processes themselves. Therefore, unless in exceptional cases, avoid using `when` blocks in processes.

### 5.3 Maximizing Outputs
Maximize outputs (reasonably) and group them only if it is relevant (e.g., `samtools` process for `convert`, `sort`, and `index`, group `bam`/`bai` but not `raw_bam` and `sorted_bam`).

### 5.4 Centralizing Publications
To easily adapt (without multiplying module versions) the publication of processes and subworkflows (not the case for workflows) specifically to each workflow, centralize publication operations in the `conf/publish.conf` file.
Note: To target a specific call of a process or subworkflow, they must be imported with a unique name. In all cases, a process (or subworkflow) cannot be called twice with the same name.

## 6. Resource Management

### 6.1 Tags and Labels
It is recommended to tag processes with labels that control at least the number of CPUs and memory allocated to tasks. Labels can be found in `modules/config/process/labels.config` and can be updated if necessary.

### 6.2 Adaptive Resource Allocation
When it is necessary to access task resources within the process (e.g., specify the number of threads or memory in a tool's execution command line), these resources should be accessed via `task.mem`, `task.ncpu`, etc.
