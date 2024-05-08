from setuptools import setup, find_namespace_packages

packages = find_namespace_packages(
  exclude=[
  ]
)

setup(
  packages=packages,
  include_package_data = True
)

