### drawConnections.smk

# produces CSV files of connections

rule drawConnections:
  input:
    in1="results/GFP20mtest/frames/{id}_seg.npy"
  output:
    outcsv="results/GFP20mtest/frames/{id}.knex.csv"
  params:
    greedy=config['minFlow']['greedy'],
    maxCost=config['minFlow']['maxCost'],
    openCost=config['minFlow']['openingCost'],
    nnDist=config['minFlow']['nearestneighborDist'],
    divScore=config['minFlow']['divisionScoreThreshold']
  shell:
    """
        python bin/fatetrack/connect_frames.py \
            --i1 {input.in1}  \
            --greedy {params.greedy} \
            --maxCost {params.maxCost} \
            --openingCost {params.openCost} \
            --nearNeighborDist {params.nnDist} \
            --divisionScoreThreshold {params.divScore} \
            --out {output.outcsv} ;
    """
