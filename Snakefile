import os
import math
#import pdb

#pdb.set_trace()

configfile: "config.yml"
singularity: "telseq.overflow-fix2.withR.sif"

mem_step_size = 3400


def get_ids_for_batch(b):
    batch_size = config["batch_size"]
    beg = batch_size * int(b)
    end = beg + batch_size
    return config["sample_ids"][beg:end]


rule telseq_results:
    params:
        input = config["input_expression"]
    output:
        "per-sample/{sample_id}.telseq.out"
    threads: 2
    resources:
        mem_mb = lambda wc, attempt:  mem_step_size * attempt
    shell:
        """
        set -uo pipefail

        tmp_dir=`mktemp -d`
        tmp_out=$tmp_dir/$(basename {output})
        tmp_in=$tmp_dir/$(basename {params.input})

        set +e
        gsutil -q -u {config[billing_project]} cp {params.input} $tmp_in
        rc=$?
        echo "download rc: $rc"
        if [[ $rc == 0 ]]; then
          read_length=$(~/samtools/samtools view -T ~/freeze11-calling/resources/ref/hs38DH.fa $tmp_in | awk '{{print length($10)}}' | head -n 1000 | sort -nr | head -n1)
          echo "sample id: {wildcards.sample_id}"
          echo "read length: $read_length"
          if [[ $read_length > 0 ]]; then
            ~/samtools/samtools view -T ~/freeze11-calling/resources/ref/hs38DH.fa -u $tmp_in | singularity exec telseq.sif telseq -k 12 -u -o $tmp_out -r $read_length /dev/stdin 
            rc=$?
          else
            rc=1
          fi
        fi

        if [[ $rc == 0 ]]; then
          cp $tmp_out {output}
          rc=$?
        fi

        rm -r $tmp_dir

        exit $rc
        """


rule combined_telseq_results_batch:
    input:
        lambda wc: [rules.telseq_results.output[0].format(sample_id=id) for id in get_ids_for_batch(wc.batch)]
    output:
       "combined-topmed-frz11/batches/b{batch}.combined.csv"
    threads: 2
    resources:
        mem_mb = lambda wc, attempt:  mem_step_size * attempt
    shell:
        """
        set -euo pipefail

        tmp_dir=`mktemp -d`
        tmp_out=$tmp_dir/$(basename {output})
       
        # There is a trailing column in telseq ouput that we remove here
        (head -n1 {input[0]} | awk -F'\t' -vOFS=, '{{NF-=1; print "sample_id",$0}}'; for f in {input}; do sample_id=$(basename $f .telseq.out); tail -n+2 $f | awk -F'\t' -vOFS=, '{{NF-=1; print "'$sample_id'",$0}}'; done) > $tmp_out
      
        mv $tmp_out {output}
        """


rule combined_telseq_results_batch_fix:
    params:
        input = lambda wc: [rules.telseq_results.output[0].format(sample_id=id) for id in get_ids_for_batch(wc.batch)]
    output:
       "combined-topmed-frz11-fix/batches/b{batch}.combined.csv"
    threads: 2
    resources:
        mem_mb = lambda wc, attempt:  mem_step_size * attempt
    shell:
        """
        set -euo pipefail

        tmp_dir=`mktemp -d`
        tmp_out=$tmp_dir/$(basename {output})
        tmp_files=""

        for f in {params.input}; do
          sample_id=$(basename $f .telseq.out)
          # There is a trailing column in telseq ouput that we remove here
          (head -n1 $f | awk -F'\t' -vOFS=, '{{NF-=1; print "sample_id",$0}}'; tail -n+2 $f | awk -F'\t' -vOFS=, '{{NF-=1; print "'$sample_id'",$0}}') > $tmp_dir/$sample_id.csv
          tmp_files="$tmp_files $tmp_dir/$sample_id.csv"
        done

        singularity exec telseq.sif ./combine-telseq-results.R $tmp_files > $tmp_out

        mv $tmp_out {output}

        rm -r $tmp_dir/
        """

rule combined_telseq_results:
    input:
        ["combined-topmed-frz11-fix/batches/b" + str(b) + ".combined.csv" for b in range(0, math.ceil(len(config["sample_ids"]) / config["batch_size"]))]
    output:
        "combined-topmed-frz11-fix/all.combined.csv"
    threads: 2
    resources:
        mem_mb = lambda wc, attempt:  mem_step_size * attempt
    shell:
        """
        set -euo pipefail

        tmp_dir=`mktemp -d`
        tmp_out=$tmp_dir/$(basename {output})

        singularity exec telseq.sif ./combine-telseq-results.R {input} > $tmp_out

        mv $tmp_out {output}
        rm -r $tmp_dir/
        """

