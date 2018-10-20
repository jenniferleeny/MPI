DIRS += sqrt3
DIRS += wire_route
# the sets of directories to do various things in
BUILDDIRS = $(DIRS:%=build-%)
CLEANDIRS = $(DIRS:%=clean-%)

all: $(BUILDDIRS)

handin: clean writeup.pdf
	@tar cvf handin.tar *
	@echo "***"
	@echo "*** SUCCESS: Created handin.tar"
	@echo "***"

writeup.pdf:
	@echo "***"
	@echo "*** ERROR: Missing writeup.pdf: please add writeup.pdf to this directory"
	@echo "***"
	@false

$(DIRS): $(BUILDDIRS)
$(BUILDDIRS):
	@$(MAKE) --no-print-directory -C $(@:build-%=%)

clean: $(CLEANDIRS)
	@rm handin.tar || true

$(CLEANDIRS): 
	@$(MAKE) --no-print-directory -C $(@:clean-%=%) clean

.PHONY: subdirs $(DIRS)
.PHONY: subdirs $(BUILDDIRS)
.PHONY: subdirs $(CLEANDIRS)
.PHONY: all clean
