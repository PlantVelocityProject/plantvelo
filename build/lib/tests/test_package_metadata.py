import ast
import os
import unittest


ROOT = os.path.dirname(os.path.dirname(__file__))


class PackageMetadataTests(unittest.TestCase):
    def read(self, relative_path):
        with open(os.path.join(ROOT, relative_path), "r", encoding="utf-8") as handle:
            return handle.read()

    def test_runtime_version_is_0_2_0(self):
        namespace = {}
        exec(self.read("plantvelo/_version.py"), namespace)
        self.assertEqual(namespace["__version__"], "0.2.0")

    def test_setup_version_is_0_2_0(self):
        tree = ast.parse(self.read("setup.py"))
        setup_call = next(
            node.value
            for node in tree.body
            if isinstance(node, ast.Expr)
            and isinstance(node.value, ast.Call)
            and getattr(node.value.func, "id", None) == "setup"
        )
        keywords = {keyword.arg: keyword.value for keyword in setup_call.keywords}
        self.assertEqual(ast.literal_eval(keywords["version"]), "0.2.0")

    def test_readme_documents_new_command_and_layers(self):
        readme = self.read("README.md")
        self.assertIn("--ir-mode", readme)
        self.assertIn("--ir-prior", readme)
        self.assertIn("--species", readme)
        self.assertIn("五字段", readme)
        self.assertIn("六字段", readme)
        self.assertIn("不可同时提供", readme)
        self.assertIn("prior-aware", readme)
        self.assertIn("velocyto-default-v1", readme)
        self.assertIn("spliced", readme)
        self.assertIn("unspliced", readme)
        self.assertIn("retained", readme)
        self.assertIn("ambiguous", readme)
        self.assertIn("U > R > S", readme)
        self.assertIn("S + U + A", readme)

    def test_readme_documents_same_schema_merge_only(self):
        readme = self.read("README.md")

        self.assertIn("三层", readme)
        self.assertIn("四层", readme)
        self.assertIn("禁止混合", readme)
        self.assertIn("velocyto_version", readme)

    def test_readme_documents_breaking_migration(self):
        readme = self.read("README.md")
        self.assertIn("0.1.x", readme)
        self.assertIn("不可直接比较", readme)
        self.assertIn("high-confidence-IR-aware", readme)

    def test_readme_does_not_document_removed_interfaces(self):
        readme = self.read("README.md")
        current_docs, migration = readme.split("## 从 0.1.x 迁移", 1)
        self.assertNotIn("--logic", current_docs)
        self.assertNotIn("--ir-flanking", current_docs)
        self.assertNotIn("PlantPermissive10X", current_docs)
        self.assertNotIn("PlantValidated10X", current_docs)
        self.assertIn("--logic", migration)
        self.assertIn("--ir-flanking", migration)


if __name__ == "__main__":
    unittest.main()
