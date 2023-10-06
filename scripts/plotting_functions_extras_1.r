library(RColorBrewer)

n.colours <- 40
duplication.colours <- c(colorRampPalette(brewer.pal(12, 'Paired'))(n.colours), 'grey30', 'grey70')
#duplication.colours <- c(rgb(0.6509804, 0.8078431, 0.8901961), rgb(0.1215686, 0.4705882, 0.7058824), rgb(0.6980392, 0.8745098, 0.5411765), rgb(0.2, 0.627451, 0.172549), rgb(0.9843137, 0.6039216, 0.6), rgb(0.8901961, 0.1019608, 0.1098039), rgb(0.9921569, 0.7490196, 0.4352941), rgb(1, 0.4980392, 0), rgb(0.7921569, 0.6980392, 0.8392157), rgb(0.4156863, 0.2392157, 0.6039216), rgb(1,0,1,0.7), rgb(0.6941176, 0.3490196, 0.1568627), rgb(0.5,0.5,0,0.8), rgb(0,0.5,0.5,0.8), rgb(0.5, 0, 0.5, 0.8), 'yellow', 'lightblue', 'pink', 'orange', 'lightgreen', 'grey30', 'grey70')
names(duplication.colours) <- paste('Dup', c(as.character(1:n.colours), '1a', '1b'), sep = '')

plot.coef <- function(this.sample, list.of.dups = NULL, diagnostics = c('FA','SS','FM','XC','BP'), start.pos = plotting.ranges$coeaexf[1], end.pos = plotting.ranges$coeaexf[2]){
	start.index <- which(compact.hmm$coeaexf[[this.sample]]$Position >= start.pos)[1]
	end.index <- tail(which(compact.hmm$coeaexf[[this.sample]]$Position <= end.pos), 1)
	plot(compact.hmm$coeaexf[[this.sample]]$Position[start.index : end.index], compact.hmm$coeaexf[[this.sample]]$Normalised_coverage[start.index : end.index], main = this.sample, ylim = c(0,12), type = 'n')
	# Add the genes to the plot
	these.gene.colours <- rep(c('grey90', 'grey95'), length.out = nrow(gene.coords$coeaexf))
	#abline(v = unlist(gene.coords$coeaexf[, c('start', 'end')]), col = rep(these.gene.colours, 2))
	rect(gene.coords$coeaexf$start, -1, gene.coords$coeaexf$end, 13, col = these.gene.colours, border = NA)
	# Add the data
	points(compact.hmm$coeaexf[[this.sample]]$Position[start.index : end.index], compact.hmm$coeaexf[[this.sample]]$Normalised_coverage[start.index : end.index])
	# We highlight along the x axis where the duplications are as predicted by the FA and SS reads
	if (!is.null(list.of.dups)){
		# Remove Dup0
		list.of.dups <- list.of.dups[-1]
		for (d in names(list.of.dups)[list.of.dups]){
			rect(known.cnvs$coeaexf[[d]]$BP$pos[1], -0.6, known.cnvs$coeaexf[[d]]$BP$pos[2], 0, col = duplication.colours[d], border = col)
			text(known.cnvs$coeaexf[[d]]$BP$pos[2], -0.4, d, cex = 1.2, adj = c(1,0))
		}
	}
	discordant.colours <- c(FA = 'blue', SS = 'cyan', FM = 'green', XC = 'red')
	y.ranges <- list(FA = c(0,3), SS = c(3,6), FM = c(6,9), XC = c(9,12))
	for (d in diagnostics){
		if (d == 'BP'){
			add.diagnostics(diagnostic.reads$coeaexf[[this.sample]][[d]]$CEP$Position, this.col = 'brown', yrange = c(0,12))
			add.diagnostics(diagnostic.reads$coeaexf[[this.sample]][[d]]$CSP$Position, this.col = 'pink', yrange = c(0,12))
		}
		else if (d == 'XC')
			add.diagnostics(diagnostic.reads$coeaexf[[this.sample]][[d]]$Position, this.col = discordant.colours[d], yrange = y.ranges[[d]])
		else
			add.diagnostics(diagnostic.reads$coeaexf[[this.sample]][[d]], this.col = discordant.colours[d], yrange = y.ranges[[d]])
	}
	these.cnv.states <- compact.hmm$coeaexf[[this.sample]]$CNV[start.index : end.index]
	lines(compact.hmm$coeaexf[[this.sample]]$Position[start.index : end.index], these.cnv.states , col = 2, lwd = 2)
	text(apply(gene.coords$coeaexf[, c('start', 'end')], 1, mean), 9, rownames(gene.coords$coeaexf), srt=90, adj = 0)
}

plot.all.coef <- function(list.of.samples, matrix.of.read.dups = read.based.cnvs$coeaexf, diagnostics = c('FA','SS','FM','XC','BP'), start.pos = plotting.ranges$coeaexf[1], end.pos = plotting.ranges$coeaexf[2]){
	x11()
	x.midpoint <- mean(c(start.pos, end.pos))
	i <- 1
	while(1){
		if (i < 1)
			i <- length(list.of.samples)
		this.sample <- list.of.samples[i]
		if (is.null(matrix.of.read.dups))
			plot.coef(this.sample, NULL, diagnostics, start.pos, end.pos)
		else
			plot.coef(this.sample, matrix.of.read.dups[this.sample,], diagnostics, start.pos, end.pos)
		x <- locator(1)$x
		if (x <= x.midpoint)
			i <- i-1
		else
			i <- i+1
	}
}

plot.coeg <- function(this.sample, list.of.dups = NULL, diagnostics = c('FA','SS','FM','XC','BP'), start.pos = plotting.ranges$coeaexg[1], end.pos = plotting.ranges$coeaexg[2]){
	start.index <- which(compact.hmm$coeaexg[[this.sample]]$Position >= start.pos)[1]
	end.index <- tail(which(compact.hmm$coeaexg[[this.sample]]$Position <= end.pos), 1)
	plot(compact.hmm$coeaexg[[this.sample]]$Position[start.index : end.index], compact.hmm$coeaexg[[this.sample]]$Normalised_coverage[start.index : end.index], main = this.sample, ylim = c(0,12), type = 'n')
	# Add the genes to the plot
	these.gene.colours <- rep(c('grey90', 'grey95'), length.out = nrow(gene.coords$coeaexg))
	#abline(v = unlist(gene.coords$coeaexg[, c('start', 'end')]), col = rep(these.gene.colours, 2))
	rect(gene.coords$coeaexg$start, -1, gene.coords$coeaexg$end, 13, col = these.gene.colours, border = NA)
	# Add the data
	points(compact.hmm$coeaexg[[this.sample]]$Position[start.index : end.index], compact.hmm$coeaexg[[this.sample]]$Normalised_coverage[start.index : end.index])
	# We highlight along the x axis where the duplications are as predicted by the FA and SS reads
	if (!is.null(list.of.dups)){
		# Remove Dup0
		list.of.dups <- list.of.dups[-1]
		for (d in names(list.of.dups)[list.of.dups]){
			if (d == 'Dup1'){
				rect(known.cnvs$coeaexg[[d]]$BP$pos[1], -0.6, 37295500, 0, col = duplication.colours[d], border = col)
				text(37295500, -0.4, d, cex = 1.2, adj = c(1,0))
			}
			else{
				rect(known.cnvs$coeaexg[[d]]$BP$pos[1], -0.6, known.cnvs$coeaexg[[d]]$BP$pos[2], 0, col = duplication.colours[d], border = col)
				text(known.cnvs$coeaexg[[d]]$BP$pos[2], -0.4, d, cex = 1.2, adj = c(1,0))
			}
		}
	}
	discordant.colours <- c(FA = 'blue', SS = 'cyan', FM = 'green', XC = 'red')
	y.ranges <- list(FA = c(0,3), SS = c(3,6), FM = c(6,9), XC = c(9,12))
	for (d in diagnostics){
		if (d == 'BP'){
			add.diagnostics(diagnostic.reads$coeaexg[[this.sample]][[d]]$CEP$Position, this.col = 'brown', yrange = c(0,12))
			add.diagnostics(diagnostic.reads$coeaexg[[this.sample]][[d]]$CSP$Position, this.col = 'pink', yrange = c(0,12))
		}
		else if (d == 'XC')
			add.diagnostics(diagnostic.reads$coeaexg[[this.sample]][[d]]$Position, this.col = discordant.colours[d], yrange = y.ranges[[d]])
		else
			add.diagnostics(diagnostic.reads$coeaexg[[this.sample]][[d]], this.col = discordant.colours[d], yrange = y.ranges[[d]])
	}
	these.cnv.states <- compact.hmm$coeaexg[[this.sample]]$CNV[start.index : end.index]
	lines(compact.hmm$coeaexg[[this.sample]]$Position[start.index : end.index], these.cnv.states , col = 2, lwd = 2)
	text(apply(gene.coords$coeaexg[, c('start', 'end')], 1, mean), 9, rownames(gene.coords$coeaexg), srt=90, adj = 0)
}

plot.all.coeg <- function(list.of.samples, matrix.of.read.dups = read.based.cnvs$coeaexg, diagnostics = c('FA','SS','FM','XC','BP'), start.pos = plotting.ranges$coeaexg[1], end.pos = plotting.ranges$coeaexg[2]){
	x11()
	x.midpoint <- mean(c(start.pos, end.pos))
	i <- 1
	while(1){
		if (i < 1)
			i <- length(list.of.samples)
		this.sample <- list.of.samples[i]
		if (is.null(matrix.of.read.dups))
			plot.coeg(this.sample, NULL, diagnostics, start.pos, end.pos)
		else
			plot.coeg(this.sample, matrix.of.read.dups[this.sample,], diagnostics, start.pos, end.pos)
		x <- locator(1)$x
		if (x <= x.midpoint)
			i <- i-1
		else
			i <- i+1
	}
}

# Make slight change to this function to accomodate new CNVs
plot.cyp6 <- function(this.sample, list.of.dups = NULL, diagnostics = c('FA','SS','FM','XC','BP'), start.pos = plotting.ranges$cyp6[1], end.pos = plotting.ranges$cyp6[2]){
	start.index <- which(compact.hmm$cyp6[[this.sample]]$Position >= start.pos)[1]
	end.index <- tail(which(compact.hmm$cyp6[[this.sample]]$Position <= end.pos), 1)
	plot(compact.hmm$cyp6[[this.sample]]$Position[start.index : end.index], compact.hmm$cyp6[[this.sample]]$Normalised_coverage[start.index : end.index], main = this.sample, ylim = c(0,12))
	# We highlight along the x axis where the duplications are as predicted by the FA and SS reads
	if (!is.null(list.of.dups)){
		Dup.order <- paste('Dup', c(as.character(c(37:20, 19, 14, 15, 13:11, 18, 10:7, 17, 6:2)), '1a', '1b'), sep = '')
		list.of.dups <- list.of.dups[Dup.order]
		for (d in names(list.of.dups)[list.of.dups]){
			if (d == 'Dup15'){
				rect(known.cnvs$cyp6[[d]]$BP$pos[1], -0.6, 28555300, 0, col = duplication.colours[d], border = col)
				text(28555300, -0.4, d, cex = 1.2, adj = c(1,0))
			}
			else if (d == 'Dup23'){
				rect(28450000, -0.6, known.cnvs$cyp6[[d]]$BP$pos[2], 0, col = duplication.colours[d], border = col)
				text(known.cnvs$cyp6[[d]]$BP$pos[2], -0.4, d, cex = 1.2, adj = c(1,0))
			}
			else if (d == 'Dup27'){
				rect(28496700, -0.6, 28499200, 0, col = duplication.colours[d], border = col)
				text(28499200, -0.4, d, cex = 1.2, adj = c(1,0))
			}
			else{
				rect(known.cnvs$cyp6[[d]]$BP$pos[1], -0.6, known.cnvs$cyp6[[d]]$BP$pos[2], 0, col = duplication.colours[d], border = col)
				text(known.cnvs$cyp6[[d]]$BP$pos[2], -0.4, d, cex = 1.2, adj = c(1,0))
			}
		}
	}
	discordant.colours <- c(FA = 'blue', SS = 'cyan', FM = 'green', XC = 'red')
	y.ranges <- list(FA = c(0,3), SS = c(3,6), FM = c(6,9), XC = c(9,12))
	for (d in diagnostics){
		if (d == 'BP'){
			add.diagnostics(diagnostic.reads$cyp6[[this.sample]][[d]]$CEP$Position, this.col = 'brown', yrange = c(0,12))
			add.diagnostics(diagnostic.reads$cyp6[[this.sample]][[d]]$CSP$Position, this.col = 'pink', yrange = c(0,12))
		}
		else if (d == 'XC')
			add.diagnostics(diagnostic.reads$cyp6[[this.sample]][[d]]$Position, this.col = discordant.colours[d], yrange = y.ranges[[d]])
		else
			add.diagnostics(diagnostic.reads$cyp6[[this.sample]][[d]], this.col = discordant.colours[d], yrange = y.ranges[[d]])
	}
	these.cnv.states <- compact.hmm$cyp6[[this.sample]]$CNV[start.index : end.index]
	lines(compact.hmm$cyp6[[this.sample]]$Position[start.index : end.index], these.cnv.states , col = 2, lwd = 2)
	# Add the genes to the plot
	these.gene.colours <- gene.colours[1:nrow(gene.coords$cyp6)]
	abline(v = unlist(gene.coords$cyp6[, c('start', 'end')]), col = rep(these.gene.colours, 2))
	text(apply(gene.coords$cyp6[, c('start', 'end')], 1, mean), 9, rownames(gene.coords$cyp6), srt=90, adj = 0, col = these.gene.colours)
}

