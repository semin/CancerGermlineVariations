#!/usr/bin/env ruby -w

require "date"
require "logger"
require "parallel"
require "pathname"
require "fileutils"

$logger         = Logger.new(STDOUT)
$logger.level   = Logger::DEBUG

def submit(cmd)
  #puts cmd;return
  out = `#{cmd}`.chomp
  if out =~ /Job\s{1}<(\d+)>/
    puts "Job #{$1} submitted."
    return Integer($1)
  else
    abort out
  end
end

$javaBin        = Pathname.new("/opt/java/jdk1.7.0_71/bin/java")
$gccBamDir1     = Pathname.new("/groups/kucherlapati/GCC/LevelII")
$homeDir        = Pathname.new("/home/sl279")
$baseDir        = $homeDir + "BiO/Research/GCC/Germline"
$newBamDir      = Pathname.new("/n/data1/hms/dbmi/park/semin/BiO/Research/Germline/bam")
$fastqDir       = $baseDir + "FASTQ"
$scriptDir      = $baseDir + "Scripts"
$vcfDir         = $baseDir + "VCF"
$panVcfDir      = $vcfDir + "PANCAN"; $panVcfDir.mkpath
$canVcfDir      = $vcfDir + "CANCERS"; $canVcfDir.mkpath
$oriQueueScript = $scriptDir + "RealignGccBams.scala"
#$oriQueueScript = $scriptDir + "RealignGccCrcBams.scala"
$installDir     = $homeDir + "BiO/Install"
$bwaBin         = Pathname.new("/home/sl279/gentoo/usr/bin/bwa")
$samtoolsBin    = Pathname.new("/home/sl279/gentoo/usr/bin/samtools")
$queueBin       = $installDir + "Sting/dist/Queue.jar"
$gatkBin        = $installDir + "GATK3/GenomeAnalysisTK.jar"
$gatkBundleDir  = $installDir + "GATK-bundle" + "2.8" + "b37"
$picardBinDir   = $installDir + "picard-tools"
$picardBin      = Pathname.new("/home/sl279/BiO/Install/picard-tools-1.133/picard.jar")
$refSeq         = $gatkBundleDir + "human_g1k_v37_decoy.fasta"
$dbsnp          = $gatkBundleDir + "dbsnp_138.b37.vcf"
$hapmap         = $gatkBundleDir + "hapmap_3.3.b37.vcf"
$g1kIndel       = $gatkBundleDir + "1000G_phase1.indels.b37.vcf"
$millsIndel     = $gatkBundleDir + "Mills_and_1000G_gold_standard.indels.b37.vcf"
$g1kOmni        = $gatkBundleDir + "1000G_omni2.5.b37.vcf"
$chrs           = (1..22).to_a << "X" << "Y" << "MT"
$chrchrs        = $chrs.map { |c| "chr#{c}" }
$tmpDir         = Pathname.new("/home/sl279/BiO/Temp"); $tmpDir.mkpath
#$cancerTypes    = %w[BLCA BRCA CRC ESCA HNSC LGG LUAD PRAD SKCM STAD THCA UCEC UVM]
$cancerTypes    = %w[CESC]


def realign_gcc_bams_using_queue
  $cancerTypes.each do |cancerType|
    next unless cancerType == "LUAD"
    newCancerBamDir = $newBamDir + cancerType; newCancerBamDir.mkpath
    newNormalStemFound  = {}
    oriBams             = if cancerType == "CRC"
                            Pathname.glob($gccBamDir1 + "*" + "*{sorted,merged}.bam").sort
                          else
                            Pathname.glob($gccBamDir1 + cancerType + "*" + "*sorted.bam").sort
                          end
    oriNormalBams       = oriBams.select { |b| b.basename.to_s =~ /TCGA\-\S{2}\-\S{4}-(\d{2})\S{1}/ && $1.to_i >= 10 }

    $logger.info "[#{cancerType}] #{oriNormalBams.size} original normal bams found"

    cancerTypeQueueScript = newCancerBamDir + $oriQueueScript.basename.sub_ext(".#{cancerType}.scala")
    FileUtils.cp($oriQueueScript, cancerTypeQueueScript)
    lsfout = cancerTypeQueueScript.sub_ext(".scala.lsfout")
    queueCmd = <<-CMD
      bsub \\
        -q park_unlimited \\
        -g /gcc/germ/#{cancerType.downcase} \\
        -o #{lsfout} \\
        #{$javaBin} -jar #{$queueBin} \\
          -S #{cancerTypeQueueScript} \\
          #{oriNormalBams.map { |b| "-i #{b}" }.join(' ')} \\
          -O #{newCancerBamDir} \\
          -R #{$refSeq} \\
          -D #{$dbsnp} \\
          -indels #{$g1kIndel} \\
          -indels #{$millsIndel} \\
          -bwa #{$bwaBin} \\
          -bt 2 \\
          -nt 2 \\
          -sg 20 \\
          -samtools #{$samtoolsBin} \\
          -tempDir #{$tmpDir} \\
          -logDir #{newCancerBamDir + ".qlog"} \\
          -bsub \\
          -jobQueue park_7d \\
          -run \\
          -startFromScratch
    CMD
    submit queueCmd
  end
end


def realign_gcc_bams
  $cancerTypes.each do |cancerType|
    next unless cancerType == "LUAD"
    newCancerDir = $newBamDir + cancerType; newCancerDir.mkpath
    newNormalStemFound  = {}
    oriBams             = if cancerType == "CRC"
                            Pathname.glob($gccBamDir2 + "*" + "*{sorted,merged}.bam").sort
                          else
                            Pathname.glob($gccBamDir1 + cancerType + "*" + "*sorted.bam").sort
                          end
    oriTumorBams        = oriBams.select { |b| b.basename.to_s =~ /TCGA\-\S{2}\-\S{4}-(\d{2})\S{1}/ && $1.to_i < 10 }
    oriNormalBams       = oriBams.select { |b| b.basename.to_s =~ /TCGA\-\S{2}\-\S{4}-(\d{2})\S{1}/ && $1.to_i >= 10 }
    queue               = cancerType == "CRC" ? "park_7d" : "long -W 720:0"

    $logger.info "[#{cancerType}] #{oriNormalBams.size} original normal bams found"

    oriNormalBams.each_with_index do |oriNormalBam, bi|
      oriNormalStem = oriNormalBam.basename.to_s
      dirRunId      = oriNormalBam.dirname.basename.to_s
      patientId     = $1 if oriNormalStem =~ /(TCGA-\S{2}-\S{4})/
      if (oriNormalStem =~ /^(TCGA-\S{2}-\S{4}-\S{3})\S+?\_(\S+?)(\_s\_\S+|)\./)
        #110525_SN590_0091_AC02W2ABXX/TCGA-AA-3681-10A-01D-1638-02_110525_SN590_0091_AC02W2ABXX_s_8.rg.sorted.bam
        #merged/TCGA-AA-3688-10A-01D_100907_HWUSI-EAS656_0025_s_4_1.merged.bam
        sampleId        = $1
        sampleRunId     = $2
        runId           = dirRunId !~ /^\d+/ ? sampleRunId : dirRunId
        newCancerBamDir = newCancerDir + dirRunId; newCancerBamDir.mkpath
        newNormalStem   = [sampleId, runId].join("___")
        if newNormalStemFound.has_key?(newNormalStem)
          $logger.warn "[#{cancerType}] #{newNormalStem} has been found already! Skipped."
          next
        else
          newNormalStemFound[newNormalStem] = true
          sai1                        = newCancerBamDir + "#{newNormalStem}.1.sai"
          sai2                        = newCancerBamDir + "#{newNormalStem}.2.sai"
          nameSortedOriNormalBamStem  = newCancerBamDir + "#{newNormalStem}.ori.nsorted"
          nameSortedOriNormalBam      = newCancerBamDir + "#{newNormalStem}.ori.nsorted.bam"
          newNormalBam                = newCancerBamDir + "#{newNormalStem}.bam"
          newNormalBai                = newCancerBamDir + "#{newNormalStem}.bam.bai"
          newNormalMateFixedBam       = newCancerBamDir + "#{newNormalStem}.mf.bam"
          newNormalMateFixedBai       = newCancerBamDir + "#{newNormalStem}.mf.bai"
          newNormalDedupedBam         = newCancerBamDir + "#{newNormalStem}.mf.dm.bam"
          newNormalDedupedBai         = newCancerBamDir + "#{newNormalStem}.mf.dm.bai"
          newNormalDedupedMetrics     = newCancerBamDir + "#{newNormalStem}.mf.dm.metrics"
          newNormalRealignedBam       = newCancerBamDir + "#{newNormalStem}.mf.dm.ra.bam"
          newNormalRealignedBai       = newCancerBamDir + "#{newNormalStem}.mf.dm.ra.bam.bai"
          newNormalRealignerInterval  = newCancerBamDir + "#{newNormalStem}.mf.dm.ra.intervals"
          newNormalRecalibratedBam    = newCancerBamDir + "#{newNormalStem}.mf.dm.ra.rc.bam"
          newNormalRecalibratedBai    = newCancerBamDir + "#{newNormalStem}.mf.dm.ra.rc.bam.bai"
          newNormalRecalibratedData   = newCancerBamDir + "#{newNormalStem}.mf.dm.ra.rc.data"

          #next if newNormalRecalibratedBam.exist?
          #warn "Cannot find #{newNormalRecalibratedBam}. Realign #{oriNormalBam}..."

          nsortCmdLsfOut = nameSortedOriNormalBam.sub_ext(".bam.lsfout")
          nsortCmdLsfOut.delete if nsortCmdLsfOut.exist?
          nameSortedOriNormalBamStem.delete if nameSortedOriNormalBamStem.exist?
          nsortCmd =<<-CMD
            bsub \\
              -g /gcc/germ/#{cancerType.downcase} \\
              -q #{queue} \\
              -o #{nsortCmdLsfOut} \\
              #{$samtoolsBin} sort -n #{oriNormalBam} #{nameSortedOriNormalBamStem}
          CMD
          nsortCmdId = submit(nsortCmd)

          sai1CmdLsfOut = sai1.sub_ext(".sai.lsfout")
          sai1CmdLsfOut.delete if sai1CmdLsfOut.exist?
          sai1.delete if sai1.exist?
          sai1Cmd =<<-CMD
            bsub \\
              -w 'done(#{nsortCmdId})' \\
              -g /gcc/germ/#{cancerType.downcase} \\
              -q #{queue} \\
              -R "rusage[mem=4000] span[hosts=1]" \\
              -M 4000000 \\
              -o #{sai1CmdLsfOut} \\
              "#{$bwaBin} aln #{$refSeq} -b1 #{nameSortedOriNormalBam} > #{sai1}"
          CMD
          sai1CmdId = submit(sai1Cmd)

          sai2CmdLsfOut = sai2.sub_ext(".sai.lsfout")
          sai2CmdLsfOut.delete if sai2CmdLsfOut.exist?
          sai2.delete if sai2.exist?
          sai2Cmd =<<-CMD
            bsub \\
              -w 'done(#{nsortCmdId})' \\
              -g /gcc/germ/#{cancerType.downcase} \\
              -q #{queue} \\
              -R "rusage[mem=4000] span[hosts=1]" \\
              -M 4000000 \\
              -o #{sai2CmdLsfOut} \\
              "#{$bwaBin} aln #{$refSeq} -b2 #{nameSortedOriNormalBam} > #{sai2}"
          CMD
          sai2CmdId = submit(sai2Cmd)

          bwaSampeCmdLsfOut = newNormalBam.sub_ext(".bam.lsfout")
          bwaSampeCmdLsfOut.delete if bwaSampeCmdLsfOut.exist?
          newNormalBam.delete if newNormalBam.exist?
          bwaSampeCmd =<<-CMD
            bsub \\
              -w 'done(#{sai1CmdId}) && done(#{sai2CmdId})' \\
              -g /gcc/germ/#{cancerType.downcase} \\
              -q #{queue} \\
              -R "rusage[mem=10000]" \\
              -M 10000000 \\
              -o #{bwaSampeCmdLsfOut} \\
              "#{$bwaBin} sampe #{$refSeq} #{sai1} #{sai2} #{nameSortedOriNormalBam} #{nameSortedOriNormalBam} \\
                -r \'@RG\\tID:0\\tPL:ILLUMINA\\tPU:#{runId}\\tLB:#{sampleId}\\tSM:#{sampleId}\' | \\
                #{$samtoolsBin} view -bSho #{newNormalBam} - &&
              rm -rf #{nameSortedOriNormalBam} #{sai1} #{sai2}"
          CMD
          bwaSampeCmdId = submit(bwaSampeCmd)

          picardCmdLsfOut = newNormalDedupedBam.sub_ext(".bam.lsfout")
          picardCmdLsfOut.delete if picardCmdLsfOut.exist?
          newNormalMateFixedBam.delete if newNormalMateFixedBam.exist?
          newNormalDedupedBam.delete if newNormalDedupedBam.exist?
          newNormalDedupedMetrics.delete if newNormalDedupedMetrics.exist?
          picardCmd =<<-CMD
            bsub \\
              -w 'done(#{bwaSampeCmdId})' \\
              -g /gcc/germ/#{cancerType.downcase} \\
              -q #{queue} \\
              -R "rusage[mem=15000]" \\
              -M 15000000 \\
              -o #{picardCmdLsfOut} \\
              "#{$javaBin} -Xmx10G -jar #{$picardBinDir}/FixMateInformation.jar \\
                INPUT=#{newNormalBam} \\
                OUTPUT=#{newNormalMateFixedBam} \\
                SORT_ORDER=coordinate \\
                VALIDATION_STRINGENCY=LENIENT \\
                TMP_DIR=#{$tmpDir} &&
              rm -rf #{newNormalBam} &&
              #{$javaBin} -Xmx10G -jar #{$picardBinDir}/MarkDuplicates.jar \\
                INPUT=#{newNormalMateFixedBam} \\
                OUTPUT=#{newNormalDedupedBam} \\
                METRICS_FILE=#{newNormalDedupedMetrics} \\
                ASSUME_SORTED=false \\
                CREATE_INDEX=true \\
                VALIDATION_STRINGENCY=LENIENT \\
                TMP_DIR=#{$tmpDir} &&
              rm -rf #{newNormalMateFixedBam} #{newNormalMateFixedBai}"
          CMD
          picardCmdId = submit(picardCmd)

          gatkRealingerTargetCreatorCmdLsfOut = newNormalRealignerInterval.sub_ext(".intervals.lsfout")
          gatkRealingerTargetCreatorCmdLsfOut.delete if gatkRealingerTargetCreatorCmdLsfOut.exist?
          newNormalRealignerInterval.delete if newNormalRealignerInterval.exist?
          gatkRealingerTargetCreatorCmd = <<-CMD
            bsub \\
              -w 'done(#{picardCmdId})' \\
              -g /gcc/germ/#{cancerType.downcase} \\
              -q #{queue} \\
              -R "rusage[mem=15000] span[hosts=1]" \\
              -M 15000000 \\
              -o #{gatkRealingerTargetCreatorCmdLsfOut} \\
              #{$javaBin} -Xmx10G -jar #{$gatkBin} \\
                -T RealignerTargetCreator \\
                -I #{newNormalDedupedBam} \\
                -o #{newNormalRealignerInterval} \\
                -R #{$refSeq} \\
                -fixMisencodedQuals \\
                -known #{$g1kIndel} \\
                -known #{$millsIndel}
          CMD
          gatkRealingerTargetCreatorCmdId = submit(gatkRealingerTargetCreatorCmd)

          gatkIndelRealignerCmdLsfOut = newNormalRealignedBam.sub_ext(".bam.lsfout")
          gatkIndelRealignerCmdLsfOut.delete if gatkIndelRealignerCmdLsfOut.exist?
          newNormalRealignedBam.delete if newNormalRealignedBam.exist?
          gatkIndelRealignerCmd = <<-CMD
            bsub \\
              -w 'done(#{gatkRealingerTargetCreatorCmdId})' \\
              -g /gcc/germ/#{cancerType.downcase} \\
              -q #{queue} \\
              -R "rusage[mem=15000]" \\
              -M 15000000 \\
              -o #{gatkIndelRealignerCmdLsfOut} \\
              "#{$javaBin} -Xmx10G -jar #{$gatkBin} \\
                -T IndelRealigner \\
                -targetIntervals #{newNormalRealignerInterval} \\
                -I #{newNormalDedupedBam} \\
                -o #{newNormalRealignedBam} \\
                -R #{$refSeq} \\
                -maxReads 100000 \\
                -fixMisencodedQuals \\
                -known #{$g1kIndel} \\
                -known #{$millsIndel} &&
              rm -rf #{newNormalDedupedBam} #{newNormalDedupedBai} &&
              samtools index #{newNormalRealignedBam}"
          CMD
          gatkIndelRealignerCmdId = submit(gatkIndelRealignerCmd)

          gatkBaseRecalibratorCmdLsfOut = newNormalRecalibratedData.sub_ext(".data.lsfout")
          gatkBaseRecalibratorCmdLsfOut.delete if gatkBaseRecalibratorCmdLsfOut.exist?
          newNormalRecalibratedData.delete if newNormalRecalibratedData.exist?
          gatkBaseRecalibratorCmd = <<-CMD
            bsub \\
              -w 'done(#{gatkIndelRealignerCmdId})' \\
              -g /gcc/germ/#{cancerType.downcase} \\
              -q #{queue} \\
              -R "rusage[mem=15000] span[hosts=1]" \\
              -M 15000000 \\
              -o #{gatkBaseRecalibratorCmdLsfOut} \\
              #{$javaBin} -Xmx10G -jar #{$gatkBin} \\
                -T BaseRecalibrator \\
                -I #{newNormalRealignedBam} \\
                -o #{newNormalRecalibratedData} \\
                -R #{$refSeq} \\
                -knownSites #{$dbsnp} \\
                -knownSites #{$g1kIndel} \\
                -knownSites #{$millsIndel}
          CMD
          gatkBaseRecalibratorCmdId = submit(gatkBaseRecalibratorCmd)

          gatkPrintReadsCmdLsfOut = newNormalRecalibratedBam.sub_ext(".bam.lsfout")
          #newNormalRecalibratedBam.delete if newNormalRecalibratedBam.exist?
          gatkPrintReadsCmd = <<-CMD
            bsub \\
              -w 'done(#{gatkBaseRecalibratorCmdId})' \\
              -g /gcc/germ/#{cancerType.downcase} \\
              -q #{queue} \\
              -R "rusage[mem=15000]" \\
              -M 15000000 \\
              -o #{gatkPrintReadsCmdLsfOut} \\
              "#{$javaBin} -Xmx10G -jar #{$gatkBin} \\
                -T PrintReads \\
                -I #{newNormalRealignedBam} \\
                -o #{newNormalRecalibratedBam} \\
                -R #{$refSeq} \\
                -BQSR #{newNormalRecalibratedData} &&
              rm -rf #{newNormalRealignedBam} #{newNormalRealignedBai} &&
              samtools index #{newNormalRecalibratedBam}"
          CMD
          gatkPrintReadsCmdId = submit(gatkPrintReadsCmd)
        end
      else
        $logger.error "[#{cancerType}] Cannot recognize #{oriNormalStem}!!!"
        exit
      end
    end
  end
end

def refine_gcc_bams
  panBamListFile = $panVcfDir + "BamFileListForGATK-PanCan.list"
  panBamListFile.readlines.each_with_index do |line, li|
    next if li != 0
    oriBam            = Pathname.new(line.strip)
    dedupedBam        = $newBamDir + oriBam.sub_ext(".dm.bam")
    dedupedBai        = dedupedBam.sub_ext(".bai")
    dedupedMetrics    = dedupedBam.sub_ext(".metrics")
    realignedBam      = dedupedBam.sub_ext(".ra.bam")
    realignedBai      = realignedBam.sub_ext(".bam.bai")
    realignerInterval = realignedBam.sub_ext(".interval")
    recalibratedBam   = realignedBam.sub_ext(".rc.bam")
    recalibratedData  = recalibratedBam.sub_ext(".data")

    # Removing PCR duplications
    picardCmdLsfOut  = dedupedBam.sub_ext(".bam.lsfout")
    picardCmd =<<-CMD
      bsub \\
        -g /gcc/germ/refine \\
        -q long -W 700:0 \\
        -R "rusage[mem=12000]" \\
        -M 12000000 \\
        -o #{picardCmdLsfOut} \\
        "#{$javaBin} -Xmx10G -jar #{$picardBinDir}/MarkDuplicates.jar \\
          INPUT=#{oriBam} \\
          OUTPUT=#{dedupedBam} \\
          METRICS_FILE=#{dedupedMetrics} \\
          ASSUME_SORTED=false \\
          SORT_ORDER=coordinate \\
          CREATE_INDEX=true \\
          VALIDATION_STRINGENCY=LENIENT \\
          TMP_DIR=#{$tmpDir}"
    CMD
    picardCmdId = submit(picardCmd)

    gatkRealingerTargetCreatorCmdLsfOut = realignerInterval.sub_ext(".interval.lsfout")
    gatkRealingerTargetCreatorCmd = <<-CMD
      bsub \\
        -w 'done(#{picardCmdId})' \\
        -g /gcc/germ/refine \\
        -q long -W 700:0 \\
        -R "rusage[mem=12000] span[hosts=1]" \\
        -M 12000000 \\
        -n 2 \\
        -o #{gatkRealingerTargetCreatorCmdLsfOut} \\
        #{$javaBin} -Xmx10G -jar #{$gatkBin} \\
          -T RealignerTargetCreator \\
          -I #{dedupedBam} \\
          -o #{realignerInterval} \\
          -R #{$refSeq} \\
          -fixMisencodedQuals \\
          -known #{$g1kIndel} \\
          -known #{$millsIndel}
    CMD
    gatkRealingerTargetCreatorCmdId = submit(gatkRealingerTargetCreatorCmd)

    gatkIndelRealignerCmdLsfOut = realignedBam.sub_ext(".bam.lsfout")
    gatkIndelRealignerCmd = <<-CMD
      bsub \\
        -w 'done(#{gatkRealingerTargetCreatorCmdId})' \\
        -g /gcc/germ/refine \\
        -q long -W 700:0 \\
        -R "rusage[mem=12000]" \\
        -M 12000000 \\
        -o #{gatkIndelRealignerCmdLsfOut} \\
        "#{$javaBin} -Xmx10G -jar #{$gatkBin} \\
          -T IndelRealigner \\
          -targetIntervals #{realignerInterval} \\
          -I #{dedupedBam} \\
          -o #{realignedBam} \\
          -R #{$refSeq} \\
          -maxReads 100000 \\
          -fixMisencodedQuals \\
          -known #{$g1kIndel} \\
          -known #{$millsIndel} &&
        rm -rf #{dedupedBam} #{dedupedBai} &&
        samtools index #{realignedBam}"
    CMD
    gatkIndelRealignerCmdId = submit(gatkIndelRealignerCmd)

    gatkBaseRecalibratorCmdLsfOut = recalibratedData.sub_ext(".data.lsfout")
    gatkBaseRecalibratorCmd = <<-CMD
      bsub \\
        -w 'done(#{gatkIndelRealignerCmdId})' \\
        -g /gcc/germ/refine \\
        -q long -W 700:0 \\
        -R "rusage[mem=12000] span[hosts=1]" \\
        -M 12000000 \\
        -n 2 \\
        -o #{gatkBaseRecalibratorCmdLsfOut} \\
        #{$javaBin} -Xmx10G -jar #{$gatkBin} \\
          -T BaseRecalibrator \\
          -I #{realignedBam} \\
          -o #{recalibratedData} \\
          -R #{$refSeq} \\
          -knownSites #{$dbsnp} \\
          -knownSites #{$g1kIndel} \\
          -knownSites #{$millsIndel}
    CMD
    gatkBaseRecalibratorCmdId = submit(gatkBaseRecalibratorCmd)

    gatkPrintReadsCmdLsfOut = recalibratedBam.sub_ext(".bam.lsfout")
    gatkPrintReadsCmd = <<-CMD
      bsub \\
        -w 'done(#{gatkBaseRecalibratorCmdId})' \\
        -g /gcc/germ/refine \\
        -q long -W 700:0 \\
        -R "rusage[mem=12000]" \\
        -M 12000000 \\
        -o #{gatkPrintReadsCmdLsfOut} \\
        "#{$javaBin} -Xmx10G -jar #{$gatkBin} \\
          -T PrintReads \\
          -I #{realignedBam} \\
          -o #{recalibratedBam} \\
          -R #{$refSeq} \\
          -BQSR #{recalibratedData} &&
        rm -rf #{realignedBam} #{realignedBai} &&
        samtools index #{recalibratedBam}"
    CMD
    submit(gatkPrintReadsCmd)
  end
end

def name_sort_bams
  $cancerTypes.each do |cancerType|
    newCancerDir = $newBamDir + cancerType; newCancerDir.mkpath
    newNormalStemFound  = {}
    oriBams = if cancerType == "CRC"
                Pathname.glob($gccBamDir1 + "*" + "*recal.bam").sort
              else
                Pathname.glob($gccBamDir1 + cancerType + "*" + "*sorted.bam").sort
              end
    oriNormalBams = oriBams.select { |b| b.basename.to_s =~ /TCGA\-\S{2}\-\S{4}-(\d{2})\S{1}/ && $1.to_i >= 10 }
    $logger.info "[#{cancerType}] #{oriNormalBams.size} original normal bams found"

    oriNormalBams.each do |oriNormalBam|
      nameSortedOriNormalBamStem = newCancerDir + oriNormalBam.basename.sub_ext(".nsorted")
      nameSortedOriNormalBam = newCancerDir + oriNormalBam.basename.sub_ext(".nsorted.bam")
      nsortCmdLsfOut = nameSortedOriNormalBam.sub_ext(".bam.lsfout")
      next if nsortCmdLsfOut.exist?
      cmd =<<-CMD
        bsub \\
          -g /germ/sort \\
          -q short -W 12:0 \\
          -o #{nsortCmdLsfOut} \\
          #{$samtoolsBin} sort -n -o #{nameSortedOriNormalBam} -T #{nameSortedOriNormalBamStem} #{oriNormalBam}
      CMD
      submit cmd
    end
  end
end

def run_bwa_mem
  ncores = 4
  $cancerTypes.each do |cancerType|
    oriBams = if cancerType == "CRC"
                Pathname.glob($gccBamDir1 + "*" + "*recal.bam").sort
              else
                Pathname.glob($gccBamDir1 + cancerType + "*" + "*sorted.bam").sort
              end
    oriNormalBams = oriBams.select { |b| b.basename.to_s =~ /TCGA\-\S{2}\-\S{4}-(\d{2})\S{1}/ && $1.to_i >= 10 }
    $logger.info "[#{cancerType}] #{oriNormalBams.size} original normal bams found"
    newCancerDir = $newBamDir + cancerType; newCancerDir.mkpath
    oriNormalBams.each do |oriNormalBam|
      oriNormalStem = oriNormalBam.basename.to_s
      dirRunId = oriNormalBam.dirname.basename.to_s
      patientId = $1 if oriNormalStem =~ /(TCGA-\S{2}-\S{4})/
      sampleId, sampleRunId = $1, $2 if (oriNormalStem =~ /^(TCGA-\S{2}-\S{4}-\S{3})\S+?\_(\S+?)(\_s\_\S+|)\./)
      #110525_SN590_0091_AC02W2ABXX/TCGA-AA-3681-10A-01D-1638-02_110525_SN590_0091_AC02W2ABXX_s_8.rg.sorted.bam
      #merged/TCGA-AA-3688-10A-01D_100907_HWUSI-EAS656_0025_s_4_1.merged.bam
      runId = dirRunId !~ /^\d+/ ? sampleRunId : dirRunId
      nameSortedOriNormalBamStem = newCancerDir + oriNormalBam.basename(".bam").sub_ext(".nsorted")
      sortedBam = nameSortedOriNormalBamStem.sub_ext(".mem.st.bam")
      lsfout = sortedBam.sub_ext(".bam.lsfout")
      raDataLsfout = sortedBam.sub_ext(".dm.ra.rc.data.lsfout")
      next if raDataLsfout.exist?
      lsfout.delete if lsfout.exist?
          #-q i2b2_1d \\
      cmd = <<-CMD
        bsub \\
          -g /germ/bwa \\
          -q mcore -W 150:0 \\
          -R "rusage[mem=5000] span[hosts=1]" \\
          -M 5000000 \\
          -n #{ncores} \\
          -o #{lsfout} "
          #{$samtoolsBin} sort -@ #{ncores} -n -l 1 -O bam -T #{nameSortedOriNormalBamStem} #{oriNormalBam} | \\
          /home/sl279/gentoo/usr/bin/bamToFastq \\
            -i /dev/stdin \\
            -fq /dev/stdout \\
            -fq2 /dev/stdout | \\
          #{$bwaBin} mem -p -M -t #{ncores} -R \'@RG\\tID:0\\tPL:ILLUMINA\\tPU:#{runId}\\tLB:#{sampleId}\\tSM:#{sampleId}\\tCN:Harvard_GCC\' -v 1 #{$refSeq} - | \\
          #{$samtoolsBin} sort -@ #{ncores} -o #{sortedBam} -T #{sortedBam.sub_ext('')} -"
          
      CMD
      submit cmd
    end
  end
end

def remove_failed_bams
  lsfFiles = Pathname.glob($baseDir + "BAM" + "*" + "*.lsfout")
  lsfFiles.each do |lsfFile|
    if lsfFile.read.match(/Exited with/m)
      lsfFileStem = lsfFile.basename.to_s.split(".")[0]
      tmpFiles = Pathname.glob(lsfFile.dirname + "#{lsfFileStem}*")
      FileUtils.rm tmpFiles
    end
  end
end

def mark_duplicates
  bams = Pathname.glob($baseDir + "BAM" + "CESC" + "*.st.bam").sort
  bams.each_with_index do |raw_bam, si|
    dedup_bam = raw_bam.sub_ext(".dm.bam")
    dedup_bai = raw_bam.sub_ext(".dm.bai")
    dedup_metrics = dedup_bam.sub_ext(".metrics")
    lsfout = dedup_bam.sub_ext(".bam.lsfout")
    next if lsfout.exist?
    #lsfout.delete if lsfout.exist?
        #-q short -W 12:0 \\
    cmd = <<-CMD
      bsub \\
        -g /germ/md \\
        -q park_1d \\
        -R "rusage[mem=6000] span[hosts=1]" \\
        -M 6000000 \\
        -o #{lsfout} \\
        #{$javaBin} -Xmx5G -jar #{$picardBin} MarkDuplicates \\
          INPUT=#{raw_bam} \\
          OUTPUT=#{dedup_bam} \\
          METRICS_FILE=#{dedup_metrics} \\
          CREATE_INDEX=true \\
          TMP_DIR=#{$tmpDir} \\
          MAX_RECORDS_IN_RAM=150000 \\
          VALIDATION_STRINGENCY=SILENT
    CMD
    submit cmd
  end
end

def remove_duplicates
  bams = Pathname.glob($baseDir + "BAM" + "CESC" + "*.st.bam").sort
  bams.each_with_index do |raw_bam, si|
    dedup_bam = raw_bam.sub_ext(".dm.bam")
    dedup_bai = raw_bam.sub_ext(".dm.bai")
    lsfout = dedup_bam.sub_ext(".bam.lsfout")
    next if lsfout.exist?
    cmd = <<-CMD
      bsub \\
        -g /germ/md \\
        -q short -W 12:0 \\
        -R "rusage[mem=6000] span[hosts=1]" \\
        -M 6000000 \\
        -o #{lsfout} \\
        /n/data1/hms/dbmi/park/semin/BiO/Install/samtools-0.1.19/samtools rmdup #{raw_bam} #{dedup_bam}
    CMD
    puts cmd
  end
end

def create_realigner_target
  bams = Pathname.glob($baseDir + "BAM" + "*" + "*.dm.bam").sort
  bams.each_with_index do |bam, si|
    realigner_interval = bam.sub_ext(".ra.intervals")
    lsfout = realigner_interval.sub_ext(".intervals.lsfout")
    next if lsfout.exist?
    #lsfout.delete if lsfout.exist?
    cmd = <<-CMD
      bsub \\
        -g /germ/ra \\
        -q short -W 12:0 \\
        -o #{lsfout} \\
        java -Xmx5G -jar #{$gatkBin} \\
          -T RealignerTargetCreator \\
          -I #{bam} \\
          -o #{realigner_interval} \\
          -R #{$refSeq} \\
          -known #{$g1kIndel} \\
          -known #{$millsIndel} \\
          -allowPotentiallyMisencodedQuals
    CMD
    puts cmd
  end
end

def realign_indels
  bams = Pathname.glob($baseDir + "BAM" + "*" + "*.dm.bam").sort
  bams.each_with_index do |bam, si|
    realigner_interval = bam.sub_ext(".ra.intervals")
    realigned_bam = bam.sub_ext(".ra.bam")
    lsfout = realigned_bam.sub_ext(".bam.lsfout")
    next if lsfout.exist?
    #lsfout.delete if lsfout.exist?
    cmd = <<-CMD
      bsub \\
        -g /germ/ra \\
        -q i2b2_1d \\
        -o #{lsfout} \\
        java -Xmx5G -jar #{$gatkBin} \\
          -T IndelRealigner \\
          -targetIntervals #{realigner_interval} \\
          -I #{bam} \\
          -o #{realigned_bam} \\
          -R #{$refSeq} \\
          -allowPotentiallyMisencodedQuals \\
          --filter_bases_not_stored \\
          -known #{$g1kIndel} \\
          -known #{$millsIndel}
    CMD
    puts cmd
  end
end

def run_bqsr
  bams = Pathname.glob($baseDir + "BAM" + "*" + "*.ra.bam").sort
  Parallel.each(bams, :in_threads => 5) do |bam|
    rc_data = bam.sub_ext(".rc.data")
    lsfout = rc_data.sub_ext(".data.lsfout")
    next if lsfout.exist?
        #-q long -W 150:0 \\
    cmd = <<-CMD
      bsub \\
        -g /germ/rc \\
        -q i2b2_1d \\
        -o #{lsfout} \\
        java -Xmx5G -jar #{$gatkBin} \\
          -T BaseRecalibrator \\
          -I #{bam} \\
          -o #{rc_data} \\
          -R #{$refSeq} \\
          -knownSites #{$dbsnp} \\
          -knownSites #{$g1kIndel} \\
          -knownSites #{$millsIndel}
    CMD
          #-allowPotentiallyMisencodedQuals
          #-fixMisencodedQuals
    puts cmd
  end
end

def print_reads
  bams = Pathname.glob($baseDir + "BAM" + "*" + "*.ra.bam").sort
  Parallel.each(bams, :in_threads => 5) do |bam|
    rc_data = bam.sub_ext(".rc.data")
    rc_bam = bam.sub_ext(".rc.bam")
    lsfout = rc_bam.sub_ext(".bam.lsfout")
    next if lsfout.exist?
    cmd = <<-CMD
      bsub \\
        -g /germ/rc \\
        -q long -W 150:0 \\
        -o #{lsfout} \\
        java -Xmx5G -jar #{$gatkBin} \\
          -T PrintReads \\
          -I #{bam} \\
          -o #{rc_bam} \\
          -R #{$refSeq} \\
          -BQSR #{rc_data}
    CMD
    submit cmd
  end
end

#realign_gcc_bams_using_queue
#name_sort_bams

#run_bwa_mem
#mark_duplicates
remove_duplicates
#create_realigner_target
#realign_indels
#run_bqsr
#print_reads

#refine_gcc_bams
#realign_gcc_bams
#remove_failed_bams
