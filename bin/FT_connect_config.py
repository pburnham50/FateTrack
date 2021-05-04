pathToSegs = "/Users/pburnham/Documents/Fatetracking_results/H2BmC_10mintest/frames/"
rlg_out_dir = "H2BmC_10min_test/"

start_Frame = 0
end_frame=11

## initial search parameter (larger will generate more candidates)
track_frameExpCoeff = 1.75 #used to look for candidate future nuclei

## Cost generation, coefficients are use to weight certain costs higher or lower
openingCost = 15 #cost for a nucleus to appear
closingCost = 30 #cost for a nucleus to disappear
costIntCoefficient = 1
costSizeCoefficient = 3
costPositionCoefficient = 1

## Mitosis scoring parameters
mitosis_RangeMultiplier=1.75 #used to look for expected mitotic partners
mitosis_MaxScore=15 #higher scores are removed from consideration of splits
DivIntensityScoreReturn = 20.0 # Mitosis Score if intensity change is high (0 if low)
DivIntensityDiffThreshold = 0.25
DivSizeScoreReturn = 100.0 # Mitosis Score if size change is high (0 if low)
DivSizeDiffThreshold = 0.3
DivSizeRatio_Min = 0.7
DivSizeRatio_Max =  3.0
DivMoveScoreReturn = 100.0 # Mitosis Score if position change is high (0 if low)
