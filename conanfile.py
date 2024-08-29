from conan import ConanFile
from conan.tools.cmake import CMakeToolchain


class Recipe(ConanFile):
    settings = "os", "compiler", "build_type", "arch"
    generators = "CMakeDeps", "VirtualRunEnv"
    options = {"shared": [True, False], "fPIC": [True, False], "parallel": [False]}
    default_options = {"shared": True, "fPIC": True, "parallel": False}

    def configure(self):
        if self.options.shared:
            # fPIC might have been removed in config_options(), so we use rm_safe
            self.options.rm_safe("fPIC")

    def generate(self):
        tc = CMakeToolchain(self)
        tc.cache_variables["MEDCOUPLING_USE_MPI"] = self.options.parallel
        tc.generate()

    def layout(self):
        self.folders.generators = "conan"

    def requirements(self):
        self.requires("hdf5/1.10.5")
        self.requires("medfile/4.1.1")

    def build_requirements(self):
        self.test_requires("cppunit/1.15.1")
