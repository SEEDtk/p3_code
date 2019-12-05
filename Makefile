TOP_DIR = ../..
include $(TOP_DIR)/tools/Makefile.common

TARGET ?= /kb/deployment
DEPLOY_RUNTIME ?= /vol/kbase/runtime

all: bin

bin: $(BIN_PERL) $(BIN_PYTHON)

deploy: deploy-all
deploy-all: deploy-client 
deploy-client: deploy-libs deploy-scripts deploy-docs

deploy-service-scripts:
	export KB_TOP=$(TARGET); \
	export KB_RUNTIME=$(DEPLOY_RUNTIME); \
	export KB_PERL_PATH=$(TARGET)/lib ; \
	for src in $(SRC_SERVICE_PERL) ; do \
	        basefile=`basename $$src`; \
	        base=`basename $$src .pl`; \
	        echo install $$src $$base ; \
	        cp $$src $(TARGET)/plbin ; \
	        $(WRAP_PERL_SCRIPT) "$(TARGET)/plbin/$$basefile" $(TARGET)/services/$(SERVICE)/bin/$$base ; \
	done


deploy-dir:
	if [ ! -d $(SERVICE_DIR) ] ; then mkdir $(SERVICE_DIR) ; fi
	if [ ! -d $(SERVICE_DIR)/bin ] ; then mkdir $(SERVICE_DIR)/bin ; fi

deploy-docs: 

include $(TOP_DIR)/tools/Makefile.common.rules
