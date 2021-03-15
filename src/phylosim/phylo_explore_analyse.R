
out.file <- "explore.txt"

dd <- read.table(out.file, header=TRUE)

pars <- c("u","v","k.regul","k.pi","deg.rate","deg.effect")
pars <- c("k.regul","k.pi")

pdf("explore.pdf", width=10, height=5*length(pars))
layout(matrix(1:(2*length(pars)), ncol=2))

for (ff in c("size", "depth")) {
	for (pp in pars) {
		dpp <- dd[,pp]
		dff <- dd[,ff]
		dpp <- dpp[is.finite(dff)]
		dff <- dff[is.finite(dff)]
		dff <- dff[order(dpp)]
		dpp <- dpp[order(dpp)]
		plot(dpp, dff, xlab=pp, ylab=ff)
		ll <- loess(dff ~ dpp)
		lines(dpp, ll$fitted, lwd=2, col="red")
	}
}

dev.off()
