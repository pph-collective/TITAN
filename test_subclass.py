from test_class import TestClass


class TestSubClass(TestClass):

    sub_class_var = 1

    def __init__(self, val, val2):
        super().__init__(val)

        self.other_val = val2

    @classmethod
    def update_class_var(cls):
        super().update_class_var()

        cls.sub_class_var += 1

    def instance_update_class_var(self):
        self.sub_class_var += 1
