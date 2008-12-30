
.PHONY: .triangle

TRIANGLE_PARAMS=-peA -q30 -a0.02

all: cooling-1.1.node \
	cooling-2.1.node \
	cooling-3.1.node \
	cooling-4.1.node

%.1.node: %.poly
	# generates .1.node, .1.ele, .1.poly , and (in case -e parameter is given) .1.edge files
	triangle $(TRIANGLE_PARAMS) $^

clean:
	rm -f *.[1-9].poly *.[1-9].node *.[1-9].ele *.[1-9].edge
	rm -f *.ALUgrid

