# --- UCSF Chimera Copyright ---
# Copyright (c) 2000 Regents of the University of California.
# All rights reserved.  This software provided pursuant to a
# license agreement containing restrictions on its disclosure,
# duplication and use.  This notice must be embedded in or
# attached to all copies, including partial copies, of the
# software or any revisions or derivations thereof.
# --- UCSF Chimera Copyright ---

import chimera

class MOLE2EMO(chimera.extension.EMO):
	_icon_shown = False
	def name(self):
		return 'MOLE 2.5'
	def description(self):
		return 'locate/characterize biomacromolecular channels'
	def icon(self):
		return self.path('molelogo.png')
	def categories(self):
		if not self._icon_shown:
			self.module()
			self._icon_shown = True
		return ['Surface/Binding Analysis']
	def activate(self):
		import chimera.dialogs
		chimera.dialogs.display(
				self.module('MainWindow').MainWindow.name)
		return None

chimera.extension.manager.registerExtension(MOLE2EMO(__file__))
