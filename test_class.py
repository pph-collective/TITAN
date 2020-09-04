class TestClass:
    class_var = 0

    def __init__(self, val):
        self.val = val

    @classmethod
    def update_class_var(cls):
        cls.class_var += 1

    @staticmethod
    def static_class_method(name):
        print(f"Hello, {name}!")
