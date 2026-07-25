.PHONY: install test lint format coverage

install:
	pip uninstall -y seamm_thermochemistry 2>/dev/null; pip install -e ".[import,test]"

test:
	pytest tests/

lint:
	black --check seamm_thermochemistry tests
	flake8 seamm_thermochemistry tests

format:
	black seamm_thermochemistry tests

coverage:
	pytest --cov=seamm_thermochemistry tests/
