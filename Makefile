all: build run

build:
	docker build -f setup/dockerfile -t geos626_seis:latest .

run:
	bash setup/run_geos626_container.sh 2>&1 | tee setup/log

clean:
	rm -rf .cache .ipython .jupyter .local .bash_history .ipynb_checkpoints