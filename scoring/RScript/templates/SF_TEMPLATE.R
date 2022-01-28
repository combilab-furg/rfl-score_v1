library(randomForest)

# load the model
load('#model_path#')

# set working directory
setwd("#tmp_path#")

# input and output file name
infn = '#input_path#'
outfn = '#output_path#'

print(paste("Read input: ", infn))
# read in input as dataframe df
df = read.table(infn, header=T, stringsAsFactors = F, sep=',')

print(df)

# get features from df
feats = df[2:(#n_features# + 1)]

# predict the binding affinity
pred = round(predict(scoring_function, newdata = feats),2)

# write output
output = data.frame(pdb = df$pdb, pred = pred)

print(paste("Write input: ", outfn))
write.table(output, outfn, sep=',', row.names = F, quote = F)

print("Done")
