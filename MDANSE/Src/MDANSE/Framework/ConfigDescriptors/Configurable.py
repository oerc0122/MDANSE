from __future__ import annotations

import json
from typing import Any

from MDANSE.IO.IOUtils import MDANSEEncoder

from .AbsConfigDesc import ConfigError, Parameter


class Configurable:
    """Allows any object that derives from it to be configurable within the MDANSE framework.

    Within that framework, to be configurable, a class must:
        - Derive from this class
    """
    @property
    def configuration(self) -> dict[str, Any]:
        return {
            name: getattr(self, name)
            for name in self.descriptors
        }


    @property
    def descriptors(self) -> dict[str, Parameter]:
        return {
            name: param
            for name, param in type(self).__dict__.items()
            if isinstance(param, Parameter)
        }

    @property
    def parameters(self) -> list[str]:
        return list(self.descriptors)

    def to_json(self) -> str:
        to_json = {
            name: obj.to_json() if isinstance(obj, Parameter) else obj
            for name, obj in self.configuration.items()
        }
        from pprint import pprint
        pprint(to_json)
        return json.dumps(to_json, cls=MDANSEEncoder)

    output_configuration = to_json

    def check_status(self) -> bool:
        try:
            for param in self.parameters:
                getattr(self, param)
        except ConfigError:
            return False

        return True
