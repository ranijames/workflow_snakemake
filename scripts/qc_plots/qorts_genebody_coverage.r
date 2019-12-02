args = commandArgs(T)

input=args[1]
output=args[2]

library(QoRTs)
packageVersion("QoRTs")
decoder = read.table(input, header=T, stringsAsFactors=F)
decoder = completeAndCheckDecoder(decoder)


res = read.qc.results.data("", decoder=decoder)

byGroup.plotter = build.plotter.colorByGroup(res)
pdf(output)
makePlot.genebody(byGroup.plotter)
makePlot.legend.over("topright",byGroup.plotter);
dev.off()
