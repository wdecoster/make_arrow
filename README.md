# make_arrow

A Rust tool to create an arrow file, containing information about read length, percent identity, and mapping quality from a cram/bam file. The files created by this tool are compatible with [NanoPlot](https://github.com/wdecoster/nanoplot) and [NanoComp](https://github.com/wdecoster/nanocomp/).

## INSTALLATION

Download a ready-to-use binary for your system from the [releases](https://github.com/wdecoster/make_arrow/releases) and add it to a directory on your $PATH .  
You may have to change the file permissions to execute it with `chmod +x make_arrow`

## USAGE

```text
make_arrow [OPTIONS] <INPUT>

Arguments:
  <INPUT>  cram or bam file (or '-' for stdin)

Options:
  -t, --threads <THREADS>  Number of parallel decompression threads to use [default: 4]
  -o, --output <OUTPUT>    Output file name [default: read_metrics.arrow]
  -h, --help               Print help information
  -V, --version            Print version information
```

## CITATION

If you use this tool, please consider citing our [publication](https://academic.oup.com/bioinformatics/article/39/5/btad311/7160911).
