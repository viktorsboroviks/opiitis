.PHONY: \
	format \
	format-cpp \
	lint \
	lint-cpp \
	clean \
	distclean

all: examples

examples:

aviize:
	git clone git@github.com:viktorsboroviks/aviize.git
	cd aviize; git checkout v2.5

iestaade:
	git clone git@github.com:viktorsboroviks/iestaade.git
	cd iestaade; git checkout v4.0

opiitis:
	git clone git@github.com:viktorsboroviks/opiitis.git
	cd opiitis; git checkout v0.1

format: format-cpp

format-cpp: include/opiitis.hpp
	clang-format -i $^

lint: lint-cpp

lint-cpp: \
		include/opiitis.hpp
	cppcheck \
		--enable=warning,portability,performance \
		--enable=style,information \
		--inconclusive \
		--library=std,posix,gnu \
		--platform=unix64 \
		--language=c++ \
		--std=c++20 \
		--inline-suppr \
		--check-level=exhaustive \
		--suppress=checkersReport \
		--checkers-report=cppcheck_report.txt \
		-I./include \
		$^

clean:
	rm -rfv *.out
	rm -rfv *.txt

distclean: clean
	rm -rfv aviize
	rm -rfv iestaade
