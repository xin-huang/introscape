# Copyright 2025 Xin Huang
#
# GNU General Public License v3.0
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, please see
#
#    https://www.gnu.org/licenses/gpl-3.0.en.html


rule validate_U50:
    input:
        outliers = "results/sai/2src/1KG/nea_den/w_{w}_y_{y}_z_{z}/{allele_type}/1KG.nea_den.{approach}.w_{w}_y_{y}_z_{z}.{allele_type}.scores.U.Q.0.999.overlapping.outliers.snps.tsv",
        vcfs = expand("results/polarized_data/1KG/2src/1KG.nea_den.{approach}.chr{i}.vcf.gz", i=np.arange(1,23), allow_missing=True),
        ref = "results/processed_data/1KG/samples/1KG.ref.samples.txt",
        tgt = "results/processed_data/1KG/samples/1KG.tgt.samples.txt",
        src = "results/processed_data/1KG/samples/1KG.nea_den.samples.txt",
    output:
        txt = "results/validation/2src/1KG/nea_den/w_{w}_y_{y}_z_{z}/U50/{allele_type}/1KG.nea_den.{approach}.w_{w}_y_{y}_z_{z}.U50.{allele_type}.outliers.overlap.tsv",
    params:
        data_dir = lambda wildcards: "polarized_data" if wildcards.allele_type == "derived" else "processed_data",
        output_dir = "results/validation/2src/1KG/nea_den/w_{w}_y_{y}_z_{z}/U50/{allele_type}",
        condition = lambda wildcards: f"'($5<{wildcards.w}&&$6>0.5&&$7=={wildcards.y[1]}&&$8=={wildcards.z[1]})'" if wildcards.allele_type == "derived" else f"'($5<{wildcards.w}&&$6>0.5&&$7=={wildcards.y[1]}&&$8=={wildcards.z[1]})||($5>1-{wildcards.w}&&$6<1-0.5&&$7==1-{wildcards.y[1]}&&$8==1-{wildcards.z[1]})'",
    shell:
        """
        ref_file={input.ref}
        tgt_file={input.tgt}
        src_file={input.src}

        sed '1d' {input.outliers} |\
            awk '{{print $1"\\t"$1":"$2"-"$3}}' |\
            while read chrom region; do
                for sample_file in ${{ref_file}} ${{tgt_file}}; do
                    prefix=$(basename ${{sample_file}} | sed -E 's/.*\\.([^.]+)\\.samples\\.txt/\\1/')
                    bcftools view results/{params.data_dir}/1KG/2src/1KG.nea_den.{wildcards.approach}.chr${{chrom}}.vcf.gz \
                        -r ${{region}} -S <(awk '{{print $2}}' ${{sample_file}}) |\
                    bcftools query -f "%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/AC\\t%INFO/AN\\n" > {params.output_dir}/1KG.{wildcards.approach}_${{prefix}}_${{region}}.txt
                done
                bcftools view results/{params.data_dir}/1KG/2src/1KG.nea_den.{wildcards.approach}.chr${{chrom}}.vcf.gz \
                    -r ${{region}} -S <(head -1 ${{src_file}} | awk '{{print $2}}') |\
                    bcftools query -f "%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/AC\\t%INFO/AN\\n" > {params.output_dir}/1KG.{wildcards.approach}_nea_${{region}}.txt
                bcftools view results/{params.data_dir}/1KG/2src/1KG.nea_den.{wildcards.approach}.chr${{chrom}}.vcf.gz \
                    -r ${{region}} -S <(tail -1 ${{src_file}} | awk '{{print $2}}') |\
                    bcftools query -f "%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/AC\\t%INFO/AN\\n" > {params.output_dir}/1KG.{wildcards.approach}_den_${{region}}.txt
            done

        sed '1d' {input.outliers} |\
            awk '{{print $1"\\t"$1":"$2"-"$3}}' |\
            while read chrom region; do
                paste {params.output_dir}/1KG.{wildcards.approach}_ref_${{region}}.txt \
                      {params.output_dir}/1KG.{wildcards.approach}_tgt_${{region}}.txt \
                      {params.output_dir}/1KG.{wildcards.approach}_nea_${{region}}.txt \
                      {params.output_dir}/1KG.{wildcards.approach}_den_${{region}}.txt |\
                      awk '{{print $1, $2, $3, $4, $5/$6, $11/$12, $17/$18, $23/$24}}' OFS="\\t" > {params.output_dir}/1KG.{wildcards.approach}_${{region}}.txt
                awk {params.condition} {params.output_dir}/1KG.{wildcards.approach}_${{region}}.txt | wc -l | awk -v region=${{region}} '{{print region"\\t"$1}}'
            done | sed 's/:/\\t/' | sed 's/-/\\t/' > {output.txt}
        """


rule validate_Q95:
    input:
        outliers = "results/sai/2src/1KG/nea_den/w_{w}_y_{y}_z_{z}/{allele_type}/1KG.nea_den.{approach}.w_{w}_y_{y}_z_{z}.{allele_type}.scores.U.Q.0.999.overlapping.outliers.snps.tsv",
        vcfs = expand("results/polarized_data/1KG/2src/1KG.nea_den.{approach}.chr{i}.vcf.gz", i=np.arange(1,23), allow_missing=True),
        ref = "results/processed_data/1KG/samples/1KG.ref.samples.txt",
        tgt = "results/processed_data/1KG/samples/1KG.tgt.samples.txt",
        src = "results/processed_data/1KG/samples/1KG.nea_den.samples.txt",
    output:
        txt = "results/validation/2src/1KG/nea_den/w_{w}_y_{y}_z_{z}/Q95/{allele_type}/1KG.nea_den.{approach}.w_{w}_y_{y}_z_{z}.Q95.{allele_type}.outliers.overlap.tsv",
    params:
        data_dir = lambda wildcards: "polarized_data" if wildcards.allele_type == "derived" else "processed_data",
        output_dir = "results/validation/2src/1KG/nea_den/w_{w}_y_{y}_z_{z}/Q95/{allele_type}",
        condition = lambda wildcards: f"'($5<{wildcards.w}&&$7=={wildcards.y[1]}&&$8=={wildcards.z[1]})'" if wildcards.allele_type == "derived" else f"'($5<{wildcards.w}&&$7=={wildcards.y[1]}&&$8=={wildcards.z[1]})||($5>1-{wildcards.w}&&$7==1-{wildcards.y[1]}&&$8==1-{wildcards.z[1]})'",
    shell:
        """
        ref_file={input.ref}
        tgt_file={input.tgt}
        src_file={input.src}

        sed '1d' {input.outliers} |\
            awk '{{print $1"\\t"$1":"$2"-"$3}}' |\
            while read chrom region; do
                for sample_file in ${{ref_file}} ${{tgt_file}}; do
                    prefix=$(basename ${{sample_file}} | sed -E 's/.*\\.([^.]+)\\.samples\\.txt/\\1/')
                    bcftools view results/{params.data_dir}/1KG/2src/1KG.nea_den.{wildcards.approach}.chr${{chrom}}.vcf.gz \
                        -r ${{region}} -S <(awk '{{print $2}}' ${{sample_file}}) |\
                    bcftools query -f "%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/AC\\t%INFO/AN\\n" > {params.output_dir}/1KG.{wildcards.approach}_${{prefix}}_${{region}}.txt
                done
                bcftools view results/{params.data_dir}/1KG/2src/1KG.nea_den.{wildcards.approach}.chr${{chrom}}.vcf.gz \
                    -r ${{region}} -S <(head -1 ${{src_file}} | awk '{{print $2}}') |\
                    bcftools query -f "%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/AC\\t%INFO/AN\\n" > {params.output_dir}/1KG.{wildcards.approach}_nea_${{region}}.txt
                bcftools view results/{params.data_dir}/1KG/2src/1KG.nea_den.{wildcards.approach}.chr${{chrom}}.vcf.gz \
                    -r ${{region}} -S <(tail -1 ${{src_file}} | awk '{{print $2}}') |\
                    bcftools query -f "%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/AC\\t%INFO/AN\\n" > {params.output_dir}/1KG.{wildcards.approach}_den_${{region}}.txt
            done

        sed '1d' {input.outliers} |\
            awk '{{print $1"\\t"$1":"$2"-"$3}}' |\
            while read chrom region; do
                paste {params.output_dir}/1KG.{wildcards.approach}_ref_${{region}}.txt \
                      {params.output_dir}/1KG.{wildcards.approach}_tgt_${{region}}.txt \
                      {params.output_dir}/1KG.{wildcards.approach}_nea_${{region}}.txt \
                      {params.output_dir}/1KG.{wildcards.approach}_den_${{region}}.txt |\
                      awk '{{print $1, $2, $3, $4, $5/$6, $11/$12, $17/$18, $23/$24}}' OFS="\\t" |\
                      awk {params.condition} |\
                      awk -v w={wildcards.w} '{{ if ($5 > 1 - w) print 1 - $6; else print $6 }}' | sort -n | awk -v q=0.95 '
                      BEGIN {{ n = 0 }} 
                      {{ a[n++] = $1 }} 
                      END {{ 
                          if (n == 0) {{ 
                              print "No data" 
                              exit 
                          }}
                          if (n == 1) {{ 
                              print a[0] 
                              exit 
                          }}
                          pos = (n - 1) * q 
                          i = int(pos) 
                          frac = pos - i 
                          if (i + 1 < n) {{ 
                              v = a[i] * (1 - frac) + a[i + 1] * frac 
                          }} else {{ 
                              v = a[i]  # when q == 1.0 exactly 
                          }}
                          printf "%.6f\\n", v 
                      }}' > {params.output_dir}/1KG.{wildcards.approach}_${{region}}.txt
               awk -v r=${{region}} '{{print r"\t"$1}}' {params.output_dir}/1KG.{wildcards.approach}_${{region}}.txt
            done | sed 's/:/\\t/' | sed 's/-/\\t/' > {output.txt}
        """
