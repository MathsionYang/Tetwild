/*
 *  Copyright (c) 2012-2016, Bruno Levy
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *  this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *  this list of conditions and the following disclaimer in the documentation
 *  and/or other materials provided with the distribution.
 *  * Neither the name of the ALICE Project-Team nor the names of its
 *  contributors may be used to endorse or promote products derived from this
 *  software without specific prior written permission.
 * 
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 *  If you modify this software, you should include a notice giving the
 *  name of the person performing the modification, the date of modification,
 *  and the reason for such modification.
 *
 *  Contact: Bruno Levy
 *
 *     Bruno.Levy@inria.fr
 *     http://www.loria.fr/~levy
 *
 *     ALICE Project
 *     LORIA, INRIA Lorraine, 
 *     Campus Scientifique, BP 239
 *     54506 VANDOEUVRE LES NANCY CEDEX 
 *     FRANCE
 *
 */

#include <geogram_gfx/lua/lua_imgui.h>
#include <geogram_gfx/glup_viewer/glup_viewer_gui.h>
#include <geogram_gfx/ImGui_ext/imgui_ext.h>
#include <geogram/lua/lua_wrap.h>
#include <geogram/basic/string.h>
#include <map>

extern void LoadImguiBindings();
extern lua_State* lState;

namespace {
    using namespace GEO;

    std::map<std::string,wchar_t> font_awesome_table;
    void init_font_awesome_table(void);
    
    int wrapper_TextInput(lua_State* L) {
	
	if(
	    lua_gettop(L) != 2 &&
	    lua_gettop(L) != 3
	) {
	    return luaL_error(
		L, "'imgui.TextInput()' invalid number of arguments"
	    );
	}
	
	if(!lua_isstring(L,1)) {
	    return luaL_error(
		L, "'imgui.TextInput()' argument 1 should be a string"
	    );
	}
	
	if(!lua_isstring(L,2)) {
	    return luaL_error(
		L, "'imgui.TextInput()' argument 2 should be a string"
	    );
	}

	ImGuiInputTextFlags flags = 0;
	
	if(lua_gettop(L) == 3) {
	    if(!lua_isnumber(L,3)) {
		return luaL_error(
		    L, "'imgui.TextInput()' argument 3 should be a number"
		);
	    }
	    flags = ImGuiInputTextFlags(lua_tonumber(L,3));
	}
	
	const char* label  = lua_tostring(L,1);
	const char* str = lua_tostring(L,2);	
	static char buff[geo_imgui_string_length];
	strcpy(buff,str);
	bool result = ImGui::InputText(
	    label, buff, geo_imgui_string_length, flags
	);
	lua_pushboolean(L,result);
	lua_pushstring(L,buff);
	
	
	return 2;
    }


    int wrapper_Combo(lua_State* L) {
	if(lua_gettop(L) != 3) {
	    return luaL_error(
		L, "'imgui.Combo()' invalid number of arguments"
	    );
	}
	
	if(!lua_isstring(L,1)) {
	    return luaL_error(
		L, "'imgui.Combo()' argument should be a string"
	    );
	}
	
	if(!lua_isstring(L,2)) {
	    return luaL_error(
		L, "'imgui.Combo()' argument should be a string"
	    );
	}

	if(!lua_isstring(L,3)) {
	    return luaL_error(
		L, "'imgui.Combo()' argument should be a string"
	    );
	}
	
	const char* label = lua_tostring(L,1);
	const char* current_item = lua_tostring(L,2);
	const char* items = lua_tostring(L,3);		

	char* lua_items = (char*)alloca(strlen(items)+2);
	strcpy(lua_items,items);
	{
	    size_t n = strlen(lua_items);
	    lua_items[n] = ';';
	    lua_items[n+1] = '\0';
	}
	
	int lua_current_item=0;
	
	const char* prev_item = lua_items;
	int nb_items = 0;

	char* p = lua_items;
	while(*p != '\0') {
	    if(*p == ';') {
		*p = '\0';
		if(!strcmp(prev_item, current_item)) {
		    lua_current_item = nb_items;
		}
		prev_item = p+1;
		++nb_items;
	    }
	    ++p;
	}
	*p = '\0'; // Double '\0' to indicate end of item list to lua.

	bool result = ImGui::Combo(label, &lua_current_item, lua_items);

	current_item = lua_items;
	while(lua_current_item > 0) {
	    while(*current_item) {
		++current_item;
	    }
	    ++current_item;
	    --lua_current_item;
	}

	lua_pushboolean(L, result);
	lua_pushstring(L, current_item);

	return 2;
    }

    int wrapper_ColorEdit3WithPalette(
	lua_State* L	
    ) {
	if(lua_gettop(L) != 4) {
	    return luaL_error(
		L, "'imgui.ColorEdit3WithPalette()' invalid number of arguments"
	    );
	}
	
	if(!lua_isstring(L,1)) {
	    return luaL_error(
		L, "'imgui.ColorEdit3WithPalette()' argument 1 should be a string"
	    );
	}

	if(!lua_isnumber(L,2)) {
	    return luaL_error(
		L, "'imgui.ColorEdit3WithPalette()' argument 2 should be a number"
	    );
	}

	if(!lua_isnumber(L,3)) {
	    return luaL_error(
		L, "'imgui.ColorEdit3WithPalette()' argument 3 should be a number"
	    );
	}

	if(!lua_isnumber(L,4)) {
	    return luaL_error(
		L, "'imgui.ColorEdit3WithPalette()' argument 4 should be a number"
	    );
	}

	const char* label = lua_tostring(L,1);
	
	float rgb[3];
	rgb[0] = float(lua_tonumber(L,2));
	rgb[1] = float(lua_tonumber(L,3));
	rgb[2] = float(lua_tonumber(L,4));	

	bool sel = ImGui::ColorEdit3WithPalette(
	    label, rgb
	);

	lua_pushboolean(L,sel);
	lua_pushnumber(L,double(rgb[0]));
	lua_pushnumber(L,double(rgb[1]));
	lua_pushnumber(L,double(rgb[2]));

	return 4;
    }

    int wrapper_ColorEdit4WithPalette(
	lua_State* L	
    ) {
	if(lua_gettop(L) != 5) {
	    return luaL_error(
		L, "'imgui.ColorEdit3WithPalette()' invalid number of arguments"
	    );
	}
	
	if(!lua_isstring(L,1)) {
	    return luaL_error(
		L, "'imgui.ColorEdit3WithPalette()' argument 1 should be a string"
	    );
	}

	if(!lua_isnumber(L,2)) {
	    return luaL_error(
		L, "'imgui.ColorEdit3WithPalette()' argument 2 should be a number"
	    );
	}

	if(!lua_isnumber(L,3)) {
	    return luaL_error(
		L, "'imgui.ColorEdit3WithPalette()' argument 3 should be a number"
	    );
	}

	if(!lua_isnumber(L,4)) {
	    return luaL_error(
		L, "'imgui.ColorEdit3WithPalette()' argument 4 should be a number"
	    );
	}

	if(!lua_isnumber(L,5)) {
	    return luaL_error(
		L, "'imgui.ColorEdit3WithPalette()' argument 5 should be a number"
	    );
	}
	
	const char* label = lua_tostring(L,1);
	
	float rgb[4];
	rgb[0] = float(lua_tonumber(L,2));
	rgb[1] = float(lua_tonumber(L,3));
	rgb[2] = float(lua_tonumber(L,4));
	rgb[3] = float(lua_tonumber(L,5));		

	bool sel = ImGui::ColorEdit4WithPalette(
	    label, rgb
	);

	lua_pushboolean(L,sel);
	lua_pushnumber(L,double(rgb[0]));
	lua_pushnumber(L,double(rgb[1]));
	lua_pushnumber(L,double(rgb[2]));
	lua_pushnumber(L,double(rgb[3]));	

	return 5;
    }

    
    int wrapper_OpenFileDialog(
	lua_State* L
    ) {
	if(lua_gettop(L) != 4) {
	    return luaL_error(
		L, "'imgui.OpenFileDialog()' invalid number of arguments"
	    );
	}

	if(!lua_isstring(L,1)) {
	    return luaL_error(
		L, "'imgui.OpenFileDialog()' argument 1 should be a string"
	    );
	}

	if(!lua_isstring(L,2)) {
	    return luaL_error(
		L, "'imgui.OpenFileDialog()' argument 2 should be a string"
	    );
	}

	if(!lua_isstring(L,3)) {
	    return luaL_error(
		L, "'imgui.OpenFileDialog()' argument 3 should be a string"
	    );
	}

	if(!lua_isnumber(L,4)) {
	    return luaL_error(
		L, "'imgui.OpenFileDialog()' argument 4 should be a number"
	    );
	}

	const char* label      = lua_tostring(L,1);
	const char* extensions = lua_tostring(L,2);
	const char* filename   = lua_tostring(L,3);
	ImGuiExtFileDialogFlags flags =
	    ImGuiExtFileDialogFlags(lua_tonumber(L,4));

	ImGui::OpenFileDialog(label, extensions, filename, flags);
	
	return 0;
    }

    int wrapper_FileDialog(
	lua_State* L
    ) {
	if(lua_gettop(L) != 2) {
	    return luaL_error(
		L, "'imgui.FileDialog()' invalid number of arguments"
	    );
	}

	if(!lua_isstring(L,1)) {
	    return luaL_error(
		L, "'imgui.OpenFileDialog()' argument 1 should be a string"
	    );
	}

	if(!lua_isstring(L,2)) {
	    return luaL_error(
		L, "'imgui.OpenFileDialog()' argument 2 should be a string"
	    );
	}

	const char* label      = lua_tostring(L,1);
	char filename[geo_imgui_string_length];

	const char* filename_in = lua_tostring(L,2);
	if(filename_in != nullptr) {
	    if(strlen(filename_in) > geo_imgui_string_length + 1) {
		Logger::err("ImGui") << "Max file name length exceeded"
				     << std::endl;
		return false;
	    }
	    strcpy(filename, filename_in);
	} else {
	    filename[0] = '\0';
	}
	
	bool result =
	    ImGui::FileDialog(label, filename, geo_imgui_string_length);

	lua_pushboolean(L,result);
	lua_pushstring(L, result ? filename : filename_in);
	
	return 2;
    }

    int wrapper_SetNextWindowPos(lua_State* L) {
	if(
	    lua_gettop(L) != 2 &&
	    lua_gettop(L) != 3
	) {
	    return luaL_error(
		L, "'imgui.SetNextWindowPos()' invalid number of arguments"
	    );
	}

	if(!lua_isnumber(L,1)) {
	    return luaL_error(
		L, "'imgui.SetNextWindowPos()' argument 1 should be a number"
	    );
	}

	if(!lua_isnumber(L,2)) {
	    return luaL_error(
		L, "'imgui.SetNextWindowPos()' argument 2 should be a number"
	    );
	}

	ImGuiCond cond = 0;
	if(lua_gettop(L) == 3) {
	    if(!lua_isnumber(L,3)) {
		return luaL_error(
		    L,
		    "'imgui.SetNextWindowPos()' argument 3 should be a number"
		);
	    }
	    cond = ImGuiCond(lua_tonumber(L,3));
	}

	
	ImGui::SetNextWindowPos(
	    ImVec2(float(lua_tonumber(L,1)), float(lua_tonumber(L,2))),
	    cond
	);

	return 0;
    }

    int wrapper_SetNextWindowSize(lua_State* L) {
	if(
	    lua_gettop(L) != 2 &&
	    lua_gettop(L) != 3
	) {
	    return luaL_error(
		L, "'imgui.SetNextWindowSize()' invalid number of arguments"
	    );
	}

	if(!lua_isnumber(L,1)) {
	    return luaL_error(
		L, "'imgui.SetNextWindowSize()' argument 1 should be a number"
	    );
	}

	if(!lua_isnumber(L,2)) {
	    return luaL_error(
		L, "'imgui.SetNextWindowSize()' argument 2 should be a number"
	    );
	}

	ImGuiCond cond = 0;
	if(lua_gettop(L) == 3) {
	    if(!lua_isnumber(L,3)) {
		return luaL_error(
		    L,
		    "'imgui.SetNextWindowSize()' argument 3 should be a number"
		);
	    }
	    cond = ImGuiCond(lua_tonumber(L,3));
	}

	
	ImGui::SetNextWindowSize(
	    ImVec2(float(lua_tonumber(L,1)), float(lua_tonumber(L,2))),
	    cond
	);

	return 0;
    }


    int wrapper_IsItemHovered(lua_State* L) {
	if(lua_gettop(L) != 0) {
	    return luaL_error(
		L, "'imgui.IsItemHovered()' invalid number of arguments"
	    );
	}
	lua_pushboolean(L,ImGui::IsItemHovered());
	return 1;
    }

    int wrapper_Text(lua_State* L) {
	const char* str = lua_tostring(L,1);
	ImGui::Text("%s",str);
	return 0;
    }

    int wrapper_SetTooltip(lua_State* L) {
	const char* str = lua_tostring(L,1);
	ImGui::SetTooltip("%s",str);
	return 0;
    }

    int wrapper_ShowStyleEditor(lua_State* L) {
	if(lua_gettop(L) != 0) {
	    return luaL_error(
		L, "'imgui.ShowStyleEditor()' invalid number of arguments"
	    );
	}
	ImGui::ShowStyleEditor(nullptr);
	return 0;
    }

    int wrapper_PushFont(lua_State* L) {
	if(lua_gettop(L) != 1) {
	    return luaL_error(
		L, "'imgui.PushFont()' invalid number of arguments"
	    );
	}
	if(!lua_isinteger(L,1)) {
	    return luaL_error(
		L, "'imgui.PushFont()' argument is not an integer"
	    );
	}
	int idx = int(lua_tointeger(L,1));
	if(idx < 0 || idx >= ImGui::GetIO().Fonts->Fonts.size()) {
	    return luaL_error(
		L, "'imgui.PushFont()' invalid font index"
	    );
	}
	ImGui::PushFont(ImGui::GetIO().Fonts->Fonts[idx]);
	return 0;
    }

    int wrapper_font_icon(lua_State* L) {
	if(lua_gettop(L) != 1) {
	    return luaL_error(
		L, "'imgui.font_icon()' invalid number of arguments"
	    );
	}
	if(!lua_isstring(L,1)) {
	    return luaL_error(
		L, "'imgui.font_icon()' argument is not an integer"
	    );
	}
	const char* K = lua_tostring(L,1);
	auto it = font_awesome_table.find(K);
	wchar_t result[2];
	result[0] = '\0';
	result[1] = '\0';	
	if(it != font_awesome_table.end()) {
	    result[0] = it->second;
	}
	std::string result_str = String::wchar_to_UTF8(result);
	lua_pushstring(L, result_str.c_str());
	return 1;
    }
}

void init_lua_imgui(lua_State* L) {
    lState = L;
    LoadImguiBindings();
    init_font_awesome_table();

    lua_pushinteger(L, ImGuiExtFileDialogFlags_Load);
    lua_setglobal(L,"ImGuiExtFileDialogFlags_Load");

    lua_pushinteger(L, ImGuiExtFileDialogFlags_Save);
    lua_setglobal(L,"ImGuiExtFileDialogFlags_Save");

    lua_pushinteger(L, ImGuiCond_Always);
    lua_setglobal(L,"ImGuiCond_Always");

    lua_pushinteger(L, ImGuiCond_Once);
    lua_setglobal(L,"ImGuiCond_Once");

    lua_pushinteger(L, ImGuiCond_FirstUseEver);
    lua_setglobal(L,"ImGuiCond_FirstUseEver");
    
    lua_pushinteger(L, ImGuiCond_Appearing);
    lua_setglobal(L,"ImGuiCond_Appearing");


    lua_pushinteger(L, ImGuiSelectableFlags_AllowDoubleClick);
    lua_setglobal(L, "ImGuiSelectableFlags_AllowDoubleClick");    
    
    lua_getglobal(L, "imgui");

    lua_pushliteral(L,"TextInput");
    lua_pushcfunction(L,wrapper_TextInput); 
    lua_settable(L,-3);

    lua_pushliteral(L,"Combo");
    lua_pushcfunction(L,wrapper_Combo); 
    lua_settable(L,-3);

    lua_pushliteral(L,"ColorEdit3WithPalette");
    lua_pushcfunction(L,wrapper_ColorEdit3WithPalette); 
    lua_settable(L,-3);

    lua_pushliteral(L,"ColorEdit4WithPalette");
    lua_pushcfunction(L,wrapper_ColorEdit4WithPalette); 
    lua_settable(L,-3);
    
    lua_pushliteral(L,"OpenFileDialog");
    lua_pushcfunction(L,wrapper_OpenFileDialog);
    lua_settable(L,-3);

    lua_pushliteral(L,"FileDialog");
    lua_pushcfunction(L,wrapper_FileDialog);
    lua_settable(L,-3);

    lua_pushliteral(L,"SetNextWindowPos");
    lua_pushcfunction(L,wrapper_SetNextWindowPos);
    lua_settable(L,-3);

    lua_pushliteral(L,"SetNextWindowSize");
    lua_pushcfunction(L,wrapper_SetNextWindowSize);
    lua_settable(L,-3);
    
    lua_pushliteral(L,"IsItemHovered");
    lua_pushcfunction(L,wrapper_IsItemHovered);
    lua_settable(L,-3);

    lua_pushliteral(L,"Text");
    lua_pushcfunction(L,wrapper_Text);
    lua_settable(L,-3);

    lua_pushliteral(L,"SetTooltip");
    lua_pushcfunction(L,wrapper_SetTooltip);
    lua_settable(L,-3);

    lua_pushliteral(L,"ShowStyleEditor");
    lua_pushcfunction(L,wrapper_ShowStyleEditor);
    lua_settable(L,-3);

    lua_pushliteral(L,"PushFont");
    lua_pushcfunction(L,wrapper_PushFont);
    lua_settable(L,-3);

    lua_pushliteral(L,"font_icon");
    lua_pushcfunction(L,wrapper_font_icon);
    lua_settable(L,-3);
    
    lua_pop(L,1);
}

namespace {

/*
   Table below generated from metadata in FontAwesome distrib: 
   cat icons.yml | awk 'BEGIN { state = 0; } {
      if(state == 0) {
         key = $1;
         gsub(":","",key);
         state = 1;
      } else if(state == 1 && $1 == "unicode:") {
         printf("font_awesome_table[\"%s\"] = 0x%s;\n",key,$2)
         state = 0;
      }  
   }'
*/    
    void init_font_awesome_table() {
	static bool initialized = false;
	if(initialized) {
	    return;
	}
	initialized = true;
	font_awesome_table["500px"] = 0xf26e;
	font_awesome_table["accessible-icon"] = 0xf368;
	font_awesome_table["accusoft"] = 0xf369;
	font_awesome_table["address-book"] = 0xf2b9;
	font_awesome_table["address-card"] = 0xf2bb;
	font_awesome_table["adjust"] = 0xf042;
	font_awesome_table["adn"] = 0xf170;
	font_awesome_table["adversal"] = 0xf36a;
	font_awesome_table["affiliatetheme"] = 0xf36b;
	font_awesome_table["air-freshener"] = 0xf5d0;
	font_awesome_table["algolia"] = 0xf36c;
	font_awesome_table["align-center"] = 0xf037;
	font_awesome_table["align-justify"] = 0xf039;
	font_awesome_table["align-left"] = 0xf036;
	font_awesome_table["align-right"] = 0xf038;
	font_awesome_table["allergies"] = 0xf461;
	font_awesome_table["amazon"] = 0xf270;
	font_awesome_table["amazon-pay"] = 0xf42c;
	font_awesome_table["ambulance"] = 0xf0f9;
	font_awesome_table["american-sign-language-interpreting"] = 0xf2a3;
	font_awesome_table["amilia"] = 0xf36d;
	font_awesome_table["anchor"] = 0xf13d;
	font_awesome_table["android"] = 0xf17b;
	font_awesome_table["angellist"] = 0xf209;
	font_awesome_table["angle-double-down"] = 0xf103;
	font_awesome_table["angle-double-left"] = 0xf100;
	font_awesome_table["angle-double-right"] = 0xf101;
	font_awesome_table["angle-double-up"] = 0xf102;
	font_awesome_table["angle-down"] = 0xf107;
	font_awesome_table["angle-left"] = 0xf104;
	font_awesome_table["angle-right"] = 0xf105;
	font_awesome_table["angle-up"] = 0xf106;
	font_awesome_table["angry"] = 0xf556;
	font_awesome_table["angrycreative"] = 0xf36e;
	font_awesome_table["angular"] = 0xf420;
	font_awesome_table["app-store"] = 0xf36f;
	font_awesome_table["app-store-ios"] = 0xf370;
	font_awesome_table["apper"] = 0xf371;
	font_awesome_table["apple"] = 0xf179;
	font_awesome_table["apple-alt"] = 0xf5d1;
	font_awesome_table["apple-pay"] = 0xf415;
	font_awesome_table["archive"] = 0xf187;
	font_awesome_table["archway"] = 0xf557;
	font_awesome_table["arrow-alt-circle-down"] = 0xf358;
	font_awesome_table["arrow-alt-circle-left"] = 0xf359;
	font_awesome_table["arrow-alt-circle-right"] = 0xf35a;
	font_awesome_table["arrow-alt-circle-up"] = 0xf35b;
	font_awesome_table["arrow-circle-down"] = 0xf0ab;
	font_awesome_table["arrow-circle-left"] = 0xf0a8;
	font_awesome_table["arrow-circle-right"] = 0xf0a9;
	font_awesome_table["arrow-circle-up"] = 0xf0aa;
	font_awesome_table["arrow-down"] = 0xf063;
	font_awesome_table["arrow-left"] = 0xf060;
	font_awesome_table["arrow-right"] = 0xf061;
	font_awesome_table["arrow-up"] = 0xf062;
	font_awesome_table["arrows-alt"] = 0xf0b2;
	font_awesome_table["arrows-alt-h"] = 0xf337;
	font_awesome_table["arrows-alt-v"] = 0xf338;
	font_awesome_table["assistive-listening-systems"] = 0xf2a2;
	font_awesome_table["asterisk"] = 0xf069;
	font_awesome_table["asymmetrik"] = 0xf372;
	font_awesome_table["at"] = 0xf1fa;
	font_awesome_table["atlas"] = 0xf558;
	font_awesome_table["atom"] = 0xf5d2;
	font_awesome_table["audible"] = 0xf373;
	font_awesome_table["audio-description"] = 0xf29e;
	font_awesome_table["autoprefixer"] = 0xf41c;
	font_awesome_table["avianex"] = 0xf374;
	font_awesome_table["aviato"] = 0xf421;
	font_awesome_table["award"] = 0xf559;
	font_awesome_table["aws"] = 0xf375;
	font_awesome_table["backspace"] = 0xf55a;
	font_awesome_table["backward"] = 0xf04a;
	font_awesome_table["balance-scale"] = 0xf24e;
	font_awesome_table["ban"] = 0xf05e;
	font_awesome_table["band-aid"] = 0xf462;
	font_awesome_table["bandcamp"] = 0xf2d5;
	font_awesome_table["barcode"] = 0xf02a;
	font_awesome_table["bars"] = 0xf0c9;
	font_awesome_table["baseball-ball"] = 0xf433;
	font_awesome_table["basketball-ball"] = 0xf434;
	font_awesome_table["bath"] = 0xf2cd;
	font_awesome_table["battery-empty"] = 0xf244;
	font_awesome_table["battery-full"] = 0xf240;
	font_awesome_table["battery-half"] = 0xf242;
	font_awesome_table["battery-quarter"] = 0xf243;
	font_awesome_table["battery-three-quarters"] = 0xf241;
	font_awesome_table["bed"] = 0xf236;
	font_awesome_table["beer"] = 0xf0fc;
	font_awesome_table["behance"] = 0xf1b4;
	font_awesome_table["behance-square"] = 0xf1b5;
	font_awesome_table["bell"] = 0xf0f3;
	font_awesome_table["bell-slash"] = 0xf1f6;
	font_awesome_table["bezier-curve"] = 0xf55b;
	font_awesome_table["bicycle"] = 0xf206;
	font_awesome_table["bimobject"] = 0xf378;
	font_awesome_table["binoculars"] = 0xf1e5;
	font_awesome_table["birthday-cake"] = 0xf1fd;
	font_awesome_table["bitbucket"] = 0xf171;
	font_awesome_table["bitcoin"] = 0xf379;
	font_awesome_table["bity"] = 0xf37a;
	font_awesome_table["black-tie"] = 0xf27e;
	font_awesome_table["blackberry"] = 0xf37b;
	font_awesome_table["blender"] = 0xf517;
	font_awesome_table["blind"] = 0xf29d;
	font_awesome_table["blogger"] = 0xf37c;
	font_awesome_table["blogger-b"] = 0xf37d;
	font_awesome_table["bluetooth"] = 0xf293;
	font_awesome_table["bluetooth-b"] = 0xf294;
	font_awesome_table["bold"] = 0xf032;
	font_awesome_table["bolt"] = 0xf0e7;
	font_awesome_table["bomb"] = 0xf1e2;
	font_awesome_table["bone"] = 0xf5d7;
	font_awesome_table["bong"] = 0xf55c;
	font_awesome_table["book"] = 0xf02d;
	font_awesome_table["book-open"] = 0xf518;
	font_awesome_table["book-reader"] = 0xf5da;
	font_awesome_table["bookmark"] = 0xf02e;
	font_awesome_table["bowling-ball"] = 0xf436;
	font_awesome_table["box"] = 0xf466;
	font_awesome_table["box-open"] = 0xf49e;
	font_awesome_table["boxes"] = 0xf468;
	font_awesome_table["braille"] = 0xf2a1;
	font_awesome_table["brain"] = 0xf5dc;
	font_awesome_table["briefcase"] = 0xf0b1;
	font_awesome_table["briefcase-medical"] = 0xf469;
	font_awesome_table["broadcast-tower"] = 0xf519;
	font_awesome_table["broom"] = 0xf51a;
	font_awesome_table["brush"] = 0xf55d;
	font_awesome_table["btc"] = 0xf15a;
	font_awesome_table["bug"] = 0xf188;
	font_awesome_table["building"] = 0xf1ad;
	font_awesome_table["bullhorn"] = 0xf0a1;
	font_awesome_table["bullseye"] = 0xf140;
	font_awesome_table["burn"] = 0xf46a;
	font_awesome_table["buromobelexperte"] = 0xf37f;
	font_awesome_table["bus"] = 0xf207;
	font_awesome_table["bus-alt"] = 0xf55e;
	font_awesome_table["buysellads"] = 0xf20d;
	font_awesome_table["calculator"] = 0xf1ec;
	font_awesome_table["calendar"] = 0xf133;
	font_awesome_table["calendar-alt"] = 0xf073;
	font_awesome_table["calendar-check"] = 0xf274;
	font_awesome_table["calendar-minus"] = 0xf272;
	font_awesome_table["calendar-plus"] = 0xf271;
	font_awesome_table["calendar-times"] = 0xf273;
	font_awesome_table["camera"] = 0xf030;
	font_awesome_table["camera-retro"] = 0xf083;
	font_awesome_table["cannabis"] = 0xf55f;
	font_awesome_table["capsules"] = 0xf46b;
	font_awesome_table["car"] = 0xf1b9;
	font_awesome_table["car-alt"] = 0xf5de;
	font_awesome_table["car-battery"] = 0xf5df;
	font_awesome_table["car-crash"] = 0xf5e1;
	font_awesome_table["car-side"] = 0xf5e4;
	font_awesome_table["caret-down"] = 0xf0d7;
	font_awesome_table["caret-left"] = 0xf0d9;
	font_awesome_table["caret-right"] = 0xf0da;
	font_awesome_table["caret-square-down"] = 0xf150;
	font_awesome_table["caret-square-left"] = 0xf191;
	font_awesome_table["caret-square-right"] = 0xf152;
	font_awesome_table["caret-square-up"] = 0xf151;
	font_awesome_table["caret-up"] = 0xf0d8;
	font_awesome_table["cart-arrow-down"] = 0xf218;
	font_awesome_table["cart-plus"] = 0xf217;
	font_awesome_table["cc-amazon-pay"] = 0xf42d;
	font_awesome_table["cc-amex"] = 0xf1f3;
	font_awesome_table["cc-apple-pay"] = 0xf416;
	font_awesome_table["cc-diners-club"] = 0xf24c;
	font_awesome_table["cc-discover"] = 0xf1f2;
	font_awesome_table["cc-jcb"] = 0xf24b;
	font_awesome_table["cc-mastercard"] = 0xf1f1;
	font_awesome_table["cc-paypal"] = 0xf1f4;
	font_awesome_table["cc-stripe"] = 0xf1f5;
	font_awesome_table["cc-visa"] = 0xf1f0;
	font_awesome_table["centercode"] = 0xf380;
	font_awesome_table["certificate"] = 0xf0a3;
	font_awesome_table["chalkboard"] = 0xf51b;
	font_awesome_table["chalkboard-teacher"] = 0xf51c;
	font_awesome_table["charging-station"] = 0xf5e7;
	font_awesome_table["chart-area"] = 0xf1fe;
	font_awesome_table["chart-bar"] = 0xf080;
	font_awesome_table["chart-line"] = 0xf201;
	font_awesome_table["chart-pie"] = 0xf200;
	font_awesome_table["check"] = 0xf00c;
	font_awesome_table["check-circle"] = 0xf058;
	font_awesome_table["check-double"] = 0xf560;
	font_awesome_table["check-square"] = 0xf14a;
	font_awesome_table["chess"] = 0xf439;
	font_awesome_table["chess-bishop"] = 0xf43a;
	font_awesome_table["chess-board"] = 0xf43c;
	font_awesome_table["chess-king"] = 0xf43f;
	font_awesome_table["chess-knight"] = 0xf441;
	font_awesome_table["chess-pawn"] = 0xf443;
	font_awesome_table["chess-queen"] = 0xf445;
	font_awesome_table["chess-rook"] = 0xf447;
	font_awesome_table["chevron-circle-down"] = 0xf13a;
	font_awesome_table["chevron-circle-left"] = 0xf137;
	font_awesome_table["chevron-circle-right"] = 0xf138;
	font_awesome_table["chevron-circle-up"] = 0xf139;
	font_awesome_table["chevron-down"] = 0xf078;
	font_awesome_table["chevron-left"] = 0xf053;
	font_awesome_table["chevron-right"] = 0xf054;
	font_awesome_table["chevron-up"] = 0xf077;
	font_awesome_table["child"] = 0xf1ae;
	font_awesome_table["chrome"] = 0xf268;
	font_awesome_table["church"] = 0xf51d;
	font_awesome_table["circle"] = 0xf111;
	font_awesome_table["circle-notch"] = 0xf1ce;
	font_awesome_table["clipboard"] = 0xf328;
	font_awesome_table["clipboard-check"] = 0xf46c;
	font_awesome_table["clipboard-list"] = 0xf46d;
	font_awesome_table["clock"] = 0xf017;
	font_awesome_table["clone"] = 0xf24d;
	font_awesome_table["closed-captioning"] = 0xf20a;
	font_awesome_table["cloud"] = 0xf0c2;
	font_awesome_table["cloud-download-alt"] = 0xf381;
	font_awesome_table["cloud-upload-alt"] = 0xf382;
	font_awesome_table["cloudscale"] = 0xf383;
	font_awesome_table["cloudsmith"] = 0xf384;
	font_awesome_table["cloudversify"] = 0xf385;
	font_awesome_table["cocktail"] = 0xf561;
	font_awesome_table["code"] = 0xf121;
	font_awesome_table["code-branch"] = 0xf126;
	font_awesome_table["codepen"] = 0xf1cb;
	font_awesome_table["codiepie"] = 0xf284;
	font_awesome_table["coffee"] = 0xf0f4;
	font_awesome_table["cog"] = 0xf013;
	font_awesome_table["cogs"] = 0xf085;
	font_awesome_table["coins"] = 0xf51e;
	font_awesome_table["columns"] = 0xf0db;
	font_awesome_table["comment"] = 0xf075;
	font_awesome_table["comment-alt"] = 0xf27a;
	font_awesome_table["comment-dots"] = 0xf4ad;
	font_awesome_table["comment-slash"] = 0xf4b3;
	font_awesome_table["comments"] = 0xf086;
	font_awesome_table["compact-disc"] = 0xf51f;
	font_awesome_table["compass"] = 0xf14e;
	font_awesome_table["compress"] = 0xf066;
	font_awesome_table["concierge-bell"] = 0xf562;
	font_awesome_table["connectdevelop"] = 0xf20e;
	font_awesome_table["contao"] = 0xf26d;
	font_awesome_table["cookie"] = 0xf563;
	font_awesome_table["cookie-bite"] = 0xf564;
	font_awesome_table["copy"] = 0xf0c5;
	font_awesome_table["copyright"] = 0xf1f9;
	font_awesome_table["couch"] = 0xf4b8;
	font_awesome_table["cpanel"] = 0xf388;
	font_awesome_table["creative-commons"] = 0xf25e;
	font_awesome_table["creative-commons-by"] = 0xf4e7;
	font_awesome_table["creative-commons-nc"] = 0xf4e8;
	font_awesome_table["creative-commons-nc-eu"] = 0xf4e9;
	font_awesome_table["creative-commons-nc-jp"] = 0xf4ea;
	font_awesome_table["creative-commons-nd"] = 0xf4eb;
	font_awesome_table["creative-commons-pd"] = 0xf4ec;
	font_awesome_table["creative-commons-pd-alt"] = 0xf4ed;
	font_awesome_table["creative-commons-remix"] = 0xf4ee;
	font_awesome_table["creative-commons-sa"] = 0xf4ef;
	font_awesome_table["creative-commons-sampling"] = 0xf4f0;
	font_awesome_table["creative-commons-sampling-plus"] = 0xf4f1;
	font_awesome_table["creative-commons-share"] = 0xf4f2;
	font_awesome_table["credit-card"] = 0xf09d;
	font_awesome_table["crop"] = 0xf125;
	font_awesome_table["crop-alt"] = 0xf565;
	font_awesome_table["crosshairs"] = 0xf05b;
	font_awesome_table["crow"] = 0xf520;
	font_awesome_table["crown"] = 0xf521;
	font_awesome_table["css3"] = 0xf13c;
	font_awesome_table["css3-alt"] = 0xf38b;
	font_awesome_table["cube"] = 0xf1b2;
	font_awesome_table["cubes"] = 0xf1b3;
	font_awesome_table["cut"] = 0xf0c4;
	font_awesome_table["cuttlefish"] = 0xf38c;
	font_awesome_table["d-and-d"] = 0xf38d;
	font_awesome_table["dashcube"] = 0xf210;
	font_awesome_table["database"] = 0xf1c0;
	font_awesome_table["deaf"] = 0xf2a4;
	font_awesome_table["delicious"] = 0xf1a5;
	font_awesome_table["deploydog"] = 0xf38e;
	font_awesome_table["deskpro"] = 0xf38f;
	font_awesome_table["desktop"] = 0xf108;
	font_awesome_table["deviantart"] = 0xf1bd;
	font_awesome_table["diagnoses"] = 0xf470;
	font_awesome_table["dice"] = 0xf522;
	font_awesome_table["dice-five"] = 0xf523;
	font_awesome_table["dice-four"] = 0xf524;
	font_awesome_table["dice-one"] = 0xf525;
	font_awesome_table["dice-six"] = 0xf526;
	font_awesome_table["dice-three"] = 0xf527;
	font_awesome_table["dice-two"] = 0xf528;
	font_awesome_table["digg"] = 0xf1a6;
	font_awesome_table["digital-ocean"] = 0xf391;
	font_awesome_table["digital-tachograph"] = 0xf566;
	font_awesome_table["directions"] = 0xf5eb;
	font_awesome_table["discord"] = 0xf392;
	font_awesome_table["discourse"] = 0xf393;
	font_awesome_table["divide"] = 0xf529;
	font_awesome_table["dizzy"] = 0xf567;
	font_awesome_table["dna"] = 0xf471;
	font_awesome_table["dochub"] = 0xf394;
	font_awesome_table["docker"] = 0xf395;
	font_awesome_table["dollar-sign"] = 0xf155;
	font_awesome_table["dolly"] = 0xf472;
	font_awesome_table["dolly-flatbed"] = 0xf474;
	font_awesome_table["donate"] = 0xf4b9;
	font_awesome_table["door-closed"] = 0xf52a;
	font_awesome_table["door-open"] = 0xf52b;
	font_awesome_table["dot-circle"] = 0xf192;
	font_awesome_table["dove"] = 0xf4ba;
	font_awesome_table["download"] = 0xf019;
	font_awesome_table["draft2digital"] = 0xf396;
	font_awesome_table["drafting-compass"] = 0xf568;
	font_awesome_table["draw-polygon"] = 0xf5ee;
	font_awesome_table["dribbble"] = 0xf17d;
	font_awesome_table["dribbble-square"] = 0xf397;
	font_awesome_table["dropbox"] = 0xf16b;
	font_awesome_table["drum"] = 0xf569;
	font_awesome_table["drum-steelpan"] = 0xf56a;
	font_awesome_table["drupal"] = 0xf1a9;
	font_awesome_table["dumbbell"] = 0xf44b;
	font_awesome_table["dyalog"] = 0xf399;
	font_awesome_table["earlybirds"] = 0xf39a;
	font_awesome_table["ebay"] = 0xf4f4;
	font_awesome_table["edge"] = 0xf282;
	font_awesome_table["edit"] = 0xf044;
	font_awesome_table["eject"] = 0xf052;
	font_awesome_table["elementor"] = 0xf430;
	font_awesome_table["ellipsis-h"] = 0xf141;
	font_awesome_table["ellipsis-v"] = 0xf142;
	font_awesome_table["ello"] = 0xf5f1;
	font_awesome_table["ember"] = 0xf423;
	font_awesome_table["empire"] = 0xf1d1;
	font_awesome_table["envelope"] = 0xf0e0;
	font_awesome_table["envelope-open"] = 0xf2b6;
	font_awesome_table["envelope-square"] = 0xf199;
	font_awesome_table["envira"] = 0xf299;
	font_awesome_table["equals"] = 0xf52c;
	font_awesome_table["eraser"] = 0xf12d;
	font_awesome_table["erlang"] = 0xf39d;
	font_awesome_table["ethereum"] = 0xf42e;
	font_awesome_table["etsy"] = 0xf2d7;
	font_awesome_table["euro-sign"] = 0xf153;
	font_awesome_table["exchange-alt"] = 0xf362;
	font_awesome_table["exclamation"] = 0xf12a;
	font_awesome_table["exclamation-circle"] = 0xf06a;
	font_awesome_table["exclamation-triangle"] = 0xf071;
	font_awesome_table["expand"] = 0xf065;
	font_awesome_table["expand-arrows-alt"] = 0xf31e;
	font_awesome_table["expeditedssl"] = 0xf23e;
	font_awesome_table["external-link-alt"] = 0xf35d;
	font_awesome_table["external-link-square-alt"] = 0xf360;
	font_awesome_table["eye"] = 0xf06e;
	font_awesome_table["eye-dropper"] = 0xf1fb;
	font_awesome_table["eye-slash"] = 0xf070;
	font_awesome_table["facebook"] = 0xf09a;
	font_awesome_table["facebook-f"] = 0xf39e;
	font_awesome_table["facebook-messenger"] = 0xf39f;
	font_awesome_table["facebook-square"] = 0xf082;
	font_awesome_table["fast-backward"] = 0xf049;
	font_awesome_table["fast-forward"] = 0xf050;
	font_awesome_table["fax"] = 0xf1ac;
	font_awesome_table["feather"] = 0xf52d;
	font_awesome_table["feather-alt"] = 0xf56b;
	font_awesome_table["female"] = 0xf182;
	font_awesome_table["fighter-jet"] = 0xf0fb;
	font_awesome_table["file"] = 0xf15b;
	font_awesome_table["file-alt"] = 0xf15c;
	font_awesome_table["file-archive"] = 0xf1c6;
	font_awesome_table["file-audio"] = 0xf1c7;
	font_awesome_table["file-code"] = 0xf1c9;
	font_awesome_table["file-contract"] = 0xf56c;
	font_awesome_table["file-download"] = 0xf56d;
	font_awesome_table["file-excel"] = 0xf1c3;
	font_awesome_table["file-export"] = 0xf56e;
	font_awesome_table["file-image"] = 0xf1c5;
	font_awesome_table["file-import"] = 0xf56f;
	font_awesome_table["file-invoice"] = 0xf570;
	font_awesome_table["file-invoice-dollar"] = 0xf571;
	font_awesome_table["file-medical"] = 0xf477;
	font_awesome_table["file-medical-alt"] = 0xf478;
	font_awesome_table["file-pdf"] = 0xf1c1;
	font_awesome_table["file-powerpoint"] = 0xf1c4;
	font_awesome_table["file-prescription"] = 0xf572;
	font_awesome_table["file-signature"] = 0xf573;
	font_awesome_table["file-upload"] = 0xf574;
	font_awesome_table["file-video"] = 0xf1c8;
	font_awesome_table["file-word"] = 0xf1c2;
	font_awesome_table["fill"] = 0xf575;
	font_awesome_table["fill-drip"] = 0xf576;
	font_awesome_table["film"] = 0xf008;
	font_awesome_table["filter"] = 0xf0b0;
	font_awesome_table["fingerprint"] = 0xf577;
	font_awesome_table["fire"] = 0xf06d;
	font_awesome_table["fire-extinguisher"] = 0xf134;
	font_awesome_table["firefox"] = 0xf269;
	font_awesome_table["first-aid"] = 0xf479;
	font_awesome_table["first-order"] = 0xf2b0;
	font_awesome_table["first-order-alt"] = 0xf50a;
	font_awesome_table["firstdraft"] = 0xf3a1;
	font_awesome_table["fish"] = 0xf578;
	font_awesome_table["flag"] = 0xf024;
	font_awesome_table["flag-checkered"] = 0xf11e;
	font_awesome_table["flask"] = 0xf0c3;
	font_awesome_table["flickr"] = 0xf16e;
	font_awesome_table["flipboard"] = 0xf44d;
	font_awesome_table["flushed"] = 0xf579;
	font_awesome_table["fly"] = 0xf417;
	font_awesome_table["folder"] = 0xf07b;
	font_awesome_table["folder-open"] = 0xf07c;
	font_awesome_table["font"] = 0xf031;
	font_awesome_table["font-awesome"] = 0xf2b4;
	font_awesome_table["font-awesome-alt"] = 0xf35c;
	font_awesome_table["font-awesome-flag"] = 0xf425;
	font_awesome_table["font-awesome-logo-full"] = 0xf4e6;
	font_awesome_table["fonticons"] = 0xf280;
	font_awesome_table["fonticons-fi"] = 0xf3a2;
	font_awesome_table["football-ball"] = 0xf44e;
	font_awesome_table["fort-awesome"] = 0xf286;
	font_awesome_table["fort-awesome-alt"] = 0xf3a3;
	font_awesome_table["forumbee"] = 0xf211;
	font_awesome_table["forward"] = 0xf04e;
	font_awesome_table["foursquare"] = 0xf180;
	font_awesome_table["free-code-camp"] = 0xf2c5;
	font_awesome_table["freebsd"] = 0xf3a4;
	font_awesome_table["frog"] = 0xf52e;
	font_awesome_table["frown"] = 0xf119;
	font_awesome_table["frown-open"] = 0xf57a;
	font_awesome_table["fulcrum"] = 0xf50b;
	font_awesome_table["futbol"] = 0xf1e3;
	font_awesome_table["galactic-republic"] = 0xf50c;
	font_awesome_table["galactic-senate"] = 0xf50d;
	font_awesome_table["gamepad"] = 0xf11b;
	font_awesome_table["gas-pump"] = 0xf52f;
	font_awesome_table["gavel"] = 0xf0e3;
	font_awesome_table["gem"] = 0xf3a5;
	font_awesome_table["genderless"] = 0xf22d;
	font_awesome_table["get-pocket"] = 0xf265;
	font_awesome_table["gg"] = 0xf260;
	font_awesome_table["gg-circle"] = 0xf261;
	font_awesome_table["gift"] = 0xf06b;
	font_awesome_table["git"] = 0xf1d3;
	font_awesome_table["git-square"] = 0xf1d2;
	font_awesome_table["github"] = 0xf09b;
	font_awesome_table["github-alt"] = 0xf113;
	font_awesome_table["github-square"] = 0xf092;
	font_awesome_table["gitkraken"] = 0xf3a6;
	font_awesome_table["gitlab"] = 0xf296;
	font_awesome_table["gitter"] = 0xf426;
	font_awesome_table["glass-martini"] = 0xf000;
	font_awesome_table["glass-martini-alt"] = 0xf57b;
	font_awesome_table["glasses"] = 0xf530;
	font_awesome_table["glide"] = 0xf2a5;
	font_awesome_table["glide-g"] = 0xf2a6;
	font_awesome_table["globe"] = 0xf0ac;
	font_awesome_table["globe-africa"] = 0xf57c;
	font_awesome_table["globe-americas"] = 0xf57d;
	font_awesome_table["globe-asia"] = 0xf57e;
	font_awesome_table["gofore"] = 0xf3a7;
	font_awesome_table["golf-ball"] = 0xf450;
	font_awesome_table["goodreads"] = 0xf3a8;
	font_awesome_table["goodreads-g"] = 0xf3a9;
	font_awesome_table["google"] = 0xf1a0;
	font_awesome_table["google-drive"] = 0xf3aa;
	font_awesome_table["google-play"] = 0xf3ab;
	font_awesome_table["google-plus"] = 0xf2b3;
	font_awesome_table["google-plus-g"] = 0xf0d5;
	font_awesome_table["google-plus-square"] = 0xf0d4;
	font_awesome_table["google-wallet"] = 0xf1ee;
	font_awesome_table["graduation-cap"] = 0xf19d;
	font_awesome_table["gratipay"] = 0xf184;
	font_awesome_table["grav"] = 0xf2d6;
	font_awesome_table["greater-than"] = 0xf531;
	font_awesome_table["greater-than-equal"] = 0xf532;
	font_awesome_table["grimace"] = 0xf57f;
	font_awesome_table["grin"] = 0xf580;
	font_awesome_table["grin-alt"] = 0xf581;
	font_awesome_table["grin-beam"] = 0xf582;
	font_awesome_table["grin-beam-sweat"] = 0xf583;
	font_awesome_table["grin-hearts"] = 0xf584;
	font_awesome_table["grin-squint"] = 0xf585;
	font_awesome_table["grin-squint-tears"] = 0xf586;
	font_awesome_table["grin-stars"] = 0xf587;
	font_awesome_table["grin-tears"] = 0xf588;
	font_awesome_table["grin-tongue"] = 0xf589;
	font_awesome_table["grin-tongue-squint"] = 0xf58a;
	font_awesome_table["grin-tongue-wink"] = 0xf58b;
	font_awesome_table["grin-wink"] = 0xf58c;
	font_awesome_table["grip-horizontal"] = 0xf58d;
	font_awesome_table["grip-vertical"] = 0xf58e;
	font_awesome_table["gripfire"] = 0xf3ac;
	font_awesome_table["grunt"] = 0xf3ad;
	font_awesome_table["gulp"] = 0xf3ae;
	font_awesome_table["h-square"] = 0xf0fd;
	font_awesome_table["hacker-news"] = 0xf1d4;
	font_awesome_table["hacker-news-square"] = 0xf3af;
	font_awesome_table["hackerrank"] = 0xf5f7;
	font_awesome_table["hand-holding"] = 0xf4bd;
	font_awesome_table["hand-holding-heart"] = 0xf4be;
	font_awesome_table["hand-holding-usd"] = 0xf4c0;
	font_awesome_table["hand-lizard"] = 0xf258;
	font_awesome_table["hand-paper"] = 0xf256;
	font_awesome_table["hand-peace"] = 0xf25b;
	font_awesome_table["hand-point-down"] = 0xf0a7;
	font_awesome_table["hand-point-left"] = 0xf0a5;
	font_awesome_table["hand-point-right"] = 0xf0a4;
	font_awesome_table["hand-point-up"] = 0xf0a6;
	font_awesome_table["hand-pointer"] = 0xf25a;
	font_awesome_table["hand-rock"] = 0xf255;
	font_awesome_table["hand-scissors"] = 0xf257;
	font_awesome_table["hand-spock"] = 0xf259;
	font_awesome_table["hands"] = 0xf4c2;
	font_awesome_table["hands-helping"] = 0xf4c4;
	font_awesome_table["handshake"] = 0xf2b5;
	font_awesome_table["hashtag"] = 0xf292;
	font_awesome_table["hdd"] = 0xf0a0;
	font_awesome_table["heading"] = 0xf1dc;
	font_awesome_table["headphones"] = 0xf025;
	font_awesome_table["headphones-alt"] = 0xf58f;
	font_awesome_table["headset"] = 0xf590;
	font_awesome_table["heart"] = 0xf004;
	font_awesome_table["heartbeat"] = 0xf21e;
	font_awesome_table["helicopter"] = 0xf533;
	font_awesome_table["highlighter"] = 0xf591;
	font_awesome_table["hips"] = 0xf452;
	font_awesome_table["hire-a-helper"] = 0xf3b0;
	font_awesome_table["history"] = 0xf1da;
	font_awesome_table["hockey-puck"] = 0xf453;
	font_awesome_table["home"] = 0xf015;
	font_awesome_table["hooli"] = 0xf427;
	font_awesome_table["hornbill"] = 0xf592;
	font_awesome_table["hospital"] = 0xf0f8;
	font_awesome_table["hospital-alt"] = 0xf47d;
	font_awesome_table["hospital-symbol"] = 0xf47e;
	font_awesome_table["hot-tub"] = 0xf593;
	font_awesome_table["hotel"] = 0xf594;
	font_awesome_table["hotjar"] = 0xf3b1;
	font_awesome_table["hourglass"] = 0xf254;
	font_awesome_table["hourglass-end"] = 0xf253;
	font_awesome_table["hourglass-half"] = 0xf252;
	font_awesome_table["hourglass-start"] = 0xf251;
	font_awesome_table["houzz"] = 0xf27c;
	font_awesome_table["html5"] = 0xf13b;
	font_awesome_table["hubspot"] = 0xf3b2;
	font_awesome_table["i-cursor"] = 0xf246;
	font_awesome_table["id-badge"] = 0xf2c1;
	font_awesome_table["id-card"] = 0xf2c2;
	font_awesome_table["id-card-alt"] = 0xf47f;
	font_awesome_table["image"] = 0xf03e;
	font_awesome_table["images"] = 0xf302;
	font_awesome_table["imdb"] = 0xf2d8;
	font_awesome_table["inbox"] = 0xf01c;
	font_awesome_table["indent"] = 0xf03c;
	font_awesome_table["industry"] = 0xf275;
	font_awesome_table["infinity"] = 0xf534;
	font_awesome_table["info"] = 0xf129;
	font_awesome_table["info-circle"] = 0xf05a;
	font_awesome_table["instagram"] = 0xf16d;
	font_awesome_table["internet-explorer"] = 0xf26b;
	font_awesome_table["ioxhost"] = 0xf208;
	font_awesome_table["italic"] = 0xf033;
	font_awesome_table["itunes"] = 0xf3b4;
	font_awesome_table["itunes-note"] = 0xf3b5;
	font_awesome_table["java"] = 0xf4e4;
	font_awesome_table["jedi-order"] = 0xf50e;
	font_awesome_table["jenkins"] = 0xf3b6;
	font_awesome_table["joget"] = 0xf3b7;
	font_awesome_table["joint"] = 0xf595;
	font_awesome_table["joomla"] = 0xf1aa;
	font_awesome_table["js"] = 0xf3b8;
	font_awesome_table["js-square"] = 0xf3b9;
	font_awesome_table["jsfiddle"] = 0xf1cc;
	font_awesome_table["kaggle"] = 0xf5fa;
	font_awesome_table["key"] = 0xf084;
	font_awesome_table["keybase"] = 0xf4f5;
	font_awesome_table["keyboard"] = 0xf11c;
	font_awesome_table["keycdn"] = 0xf3ba;
	font_awesome_table["kickstarter"] = 0xf3bb;
	font_awesome_table["kickstarter-k"] = 0xf3bc;
	font_awesome_table["kiss"] = 0xf596;
	font_awesome_table["kiss-beam"] = 0xf597;
	font_awesome_table["kiss-wink-heart"] = 0xf598;
	font_awesome_table["kiwi-bird"] = 0xf535;
	font_awesome_table["korvue"] = 0xf42f;
	font_awesome_table["language"] = 0xf1ab;
	font_awesome_table["laptop"] = 0xf109;
	font_awesome_table["laptop-code"] = 0xf5fc;
	font_awesome_table["laravel"] = 0xf3bd;
	font_awesome_table["lastfm"] = 0xf202;
	font_awesome_table["lastfm-square"] = 0xf203;
	font_awesome_table["laugh"] = 0xf599;
	font_awesome_table["laugh-beam"] = 0xf59a;
	font_awesome_table["laugh-squint"] = 0xf59b;
	font_awesome_table["laugh-wink"] = 0xf59c;
	font_awesome_table["layer-group"] = 0xf5fd;
	font_awesome_table["leaf"] = 0xf06c;
	font_awesome_table["leanpub"] = 0xf212;
	font_awesome_table["lemon"] = 0xf094;
	font_awesome_table["less"] = 0xf41d;
	font_awesome_table["less-than"] = 0xf536;
	font_awesome_table["less-than-equal"] = 0xf537;
	font_awesome_table["level-down-alt"] = 0xf3be;
	font_awesome_table["level-up-alt"] = 0xf3bf;
	font_awesome_table["life-ring"] = 0xf1cd;
	font_awesome_table["lightbulb"] = 0xf0eb;
	font_awesome_table["line"] = 0xf3c0;
	font_awesome_table["link"] = 0xf0c1;
	font_awesome_table["linkedin"] = 0xf08c;
	font_awesome_table["linkedin-in"] = 0xf0e1;
	font_awesome_table["linode"] = 0xf2b8;
	font_awesome_table["linux"] = 0xf17c;
	font_awesome_table["lira-sign"] = 0xf195;
	font_awesome_table["list"] = 0xf03a;
	font_awesome_table["list-alt"] = 0xf022;
	font_awesome_table["list-ol"] = 0xf0cb;
	font_awesome_table["list-ul"] = 0xf0ca;
	font_awesome_table["location-arrow"] = 0xf124;
	font_awesome_table["lock"] = 0xf023;
	font_awesome_table["lock-open"] = 0xf3c1;
	font_awesome_table["long-arrow-alt-down"] = 0xf309;
	font_awesome_table["long-arrow-alt-left"] = 0xf30a;
	font_awesome_table["long-arrow-alt-right"] = 0xf30b;
	font_awesome_table["long-arrow-alt-up"] = 0xf30c;
	font_awesome_table["low-vision"] = 0xf2a8;
	font_awesome_table["luggage-cart"] = 0xf59d;
	font_awesome_table["lyft"] = 0xf3c3;
	font_awesome_table["magento"] = 0xf3c4;
	font_awesome_table["magic"] = 0xf0d0;
	font_awesome_table["magnet"] = 0xf076;
	font_awesome_table["mailchimp"] = 0xf59e;
	font_awesome_table["male"] = 0xf183;
	font_awesome_table["mandalorian"] = 0xf50f;
	font_awesome_table["map"] = 0xf279;
	font_awesome_table["map-marked"] = 0xf59f;
	font_awesome_table["map-marked-alt"] = 0xf5a0;
	font_awesome_table["map-marker"] = 0xf041;
	font_awesome_table["map-marker-alt"] = 0xf3c5;
	font_awesome_table["map-pin"] = 0xf276;
	font_awesome_table["map-signs"] = 0xf277;
	font_awesome_table["markdown"] = 0xf60f;
	font_awesome_table["marker"] = 0xf5a1;
	font_awesome_table["mars"] = 0xf222;
	font_awesome_table["mars-double"] = 0xf227;
	font_awesome_table["mars-stroke"] = 0xf229;
	font_awesome_table["mars-stroke-h"] = 0xf22b;
	font_awesome_table["mars-stroke-v"] = 0xf22a;
	font_awesome_table["mastodon"] = 0xf4f6;
	font_awesome_table["maxcdn"] = 0xf136;
	font_awesome_table["medal"] = 0xf5a2;
	font_awesome_table["medapps"] = 0xf3c6;
	font_awesome_table["medium"] = 0xf23a;
	font_awesome_table["medium-m"] = 0xf3c7;
	font_awesome_table["medkit"] = 0xf0fa;
	font_awesome_table["medrt"] = 0xf3c8;
	font_awesome_table["meetup"] = 0xf2e0;
	font_awesome_table["megaport"] = 0xf5a3;
	font_awesome_table["meh"] = 0xf11a;
	font_awesome_table["meh-blank"] = 0xf5a4;
	font_awesome_table["meh-rolling-eyes"] = 0xf5a5;
	font_awesome_table["memory"] = 0xf538;
	font_awesome_table["mercury"] = 0xf223;
	font_awesome_table["microchip"] = 0xf2db;
	font_awesome_table["microphone"] = 0xf130;
	font_awesome_table["microphone-alt"] = 0xf3c9;
	font_awesome_table["microphone-alt-slash"] = 0xf539;
	font_awesome_table["microphone-slash"] = 0xf131;
	font_awesome_table["microscope"] = 0xf610;
	font_awesome_table["microsoft"] = 0xf3ca;
	font_awesome_table["minus"] = 0xf068;
	font_awesome_table["minus-circle"] = 0xf056;
	font_awesome_table["minus-square"] = 0xf146;
	font_awesome_table["mix"] = 0xf3cb;
	font_awesome_table["mixcloud"] = 0xf289;
	font_awesome_table["mizuni"] = 0xf3cc;
	font_awesome_table["mobile"] = 0xf10b;
	font_awesome_table["mobile-alt"] = 0xf3cd;
	font_awesome_table["modx"] = 0xf285;
	font_awesome_table["monero"] = 0xf3d0;
	font_awesome_table["money-bill"] = 0xf0d6;
	font_awesome_table["money-bill-alt"] = 0xf3d1;
	font_awesome_table["money-bill-wave"] = 0xf53a;
	font_awesome_table["money-bill-wave-alt"] = 0xf53b;
	font_awesome_table["money-check"] = 0xf53c;
	font_awesome_table["money-check-alt"] = 0xf53d;
	font_awesome_table["monument"] = 0xf5a6;
	font_awesome_table["moon"] = 0xf186;
	font_awesome_table["mortar-pestle"] = 0xf5a7;
	font_awesome_table["motorcycle"] = 0xf21c;
	font_awesome_table["mouse-pointer"] = 0xf245;
	font_awesome_table["music"] = 0xf001;
	font_awesome_table["napster"] = 0xf3d2;
	font_awesome_table["neos"] = 0xf612;
	font_awesome_table["neuter"] = 0xf22c;
	font_awesome_table["newspaper"] = 0xf1ea;
	font_awesome_table["nimblr"] = 0xf5a8;
	font_awesome_table["nintendo-switch"] = 0xf418;
	font_awesome_table["node"] = 0xf419;
	font_awesome_table["node-js"] = 0xf3d3;
	font_awesome_table["not-equal"] = 0xf53e;
	font_awesome_table["notes-medical"] = 0xf481;
	font_awesome_table["npm"] = 0xf3d4;
	font_awesome_table["ns8"] = 0xf3d5;
	font_awesome_table["nutritionix"] = 0xf3d6;
	font_awesome_table["object-group"] = 0xf247;
	font_awesome_table["object-ungroup"] = 0xf248;
	font_awesome_table["odnoklassniki"] = 0xf263;
	font_awesome_table["odnoklassniki-square"] = 0xf264;
	font_awesome_table["oil-can"] = 0xf613;
	font_awesome_table["old-republic"] = 0xf510;
	font_awesome_table["opencart"] = 0xf23d;
	font_awesome_table["openid"] = 0xf19b;
	font_awesome_table["opera"] = 0xf26a;
	font_awesome_table["optin-monster"] = 0xf23c;
	font_awesome_table["osi"] = 0xf41a;
	font_awesome_table["outdent"] = 0xf03b;
	font_awesome_table["page4"] = 0xf3d7;
	font_awesome_table["pagelines"] = 0xf18c;
	font_awesome_table["paint-brush"] = 0xf1fc;
	font_awesome_table["paint-roller"] = 0xf5aa;
	font_awesome_table["palette"] = 0xf53f;
	font_awesome_table["palfed"] = 0xf3d8;
	font_awesome_table["pallet"] = 0xf482;
	font_awesome_table["paper-plane"] = 0xf1d8;
	font_awesome_table["paperclip"] = 0xf0c6;
	font_awesome_table["parachute-box"] = 0xf4cd;
	font_awesome_table["paragraph"] = 0xf1dd;
	font_awesome_table["parking"] = 0xf540;
	font_awesome_table["passport"] = 0xf5ab;
	font_awesome_table["paste"] = 0xf0ea;
	font_awesome_table["patreon"] = 0xf3d9;
	font_awesome_table["pause"] = 0xf04c;
	font_awesome_table["pause-circle"] = 0xf28b;
	font_awesome_table["paw"] = 0xf1b0;
	font_awesome_table["paypal"] = 0xf1ed;
	font_awesome_table["pen"] = 0xf304;
	font_awesome_table["pen-alt"] = 0xf305;
	font_awesome_table["pen-fancy"] = 0xf5ac;
	font_awesome_table["pen-nib"] = 0xf5ad;
	font_awesome_table["pen-square"] = 0xf14b;
	font_awesome_table["pencil-alt"] = 0xf303;
	font_awesome_table["pencil-ruler"] = 0xf5ae;
	font_awesome_table["people-carry"] = 0xf4ce;
	font_awesome_table["percent"] = 0xf295;
	font_awesome_table["percentage"] = 0xf541;
	font_awesome_table["periscope"] = 0xf3da;
	font_awesome_table["phabricator"] = 0xf3db;
	font_awesome_table["phoenix-framework"] = 0xf3dc;
	font_awesome_table["phoenix-squadron"] = 0xf511;
	font_awesome_table["phone"] = 0xf095;
	font_awesome_table["phone-slash"] = 0xf3dd;
	font_awesome_table["phone-square"] = 0xf098;
	font_awesome_table["phone-volume"] = 0xf2a0;
	font_awesome_table["php"] = 0xf457;
	font_awesome_table["pied-piper"] = 0xf2ae;
	font_awesome_table["pied-piper-alt"] = 0xf1a8;
	font_awesome_table["pied-piper-hat"] = 0xf4e5;
	font_awesome_table["pied-piper-pp"] = 0xf1a7;
	font_awesome_table["piggy-bank"] = 0xf4d3;
	font_awesome_table["pills"] = 0xf484;
	font_awesome_table["pinterest"] = 0xf0d2;
	font_awesome_table["pinterest-p"] = 0xf231;
	font_awesome_table["pinterest-square"] = 0xf0d3;
	font_awesome_table["plane"] = 0xf072;
	font_awesome_table["plane-arrival"] = 0xf5af;
	font_awesome_table["plane-departure"] = 0xf5b0;
	font_awesome_table["play"] = 0xf04b;
	font_awesome_table["play-circle"] = 0xf144;
	font_awesome_table["playstation"] = 0xf3df;
	font_awesome_table["plug"] = 0xf1e6;
	font_awesome_table["plus"] = 0xf067;
	font_awesome_table["plus-circle"] = 0xf055;
	font_awesome_table["plus-square"] = 0xf0fe;
	font_awesome_table["podcast"] = 0xf2ce;
	font_awesome_table["poo"] = 0xf2fe;
	font_awesome_table["poop"] = 0xf619;
	font_awesome_table["portrait"] = 0xf3e0;
	font_awesome_table["pound-sign"] = 0xf154;
	font_awesome_table["power-off"] = 0xf011;
	font_awesome_table["prescription"] = 0xf5b1;
	font_awesome_table["prescription-bottle"] = 0xf485;
	font_awesome_table["prescription-bottle-alt"] = 0xf486;
	font_awesome_table["print"] = 0xf02f;
	font_awesome_table["procedures"] = 0xf487;
	font_awesome_table["product-hunt"] = 0xf288;
	font_awesome_table["project-diagram"] = 0xf542;
	font_awesome_table["pushed"] = 0xf3e1;
	font_awesome_table["puzzle-piece"] = 0xf12e;
	font_awesome_table["python"] = 0xf3e2;
	font_awesome_table["qq"] = 0xf1d6;
	font_awesome_table["qrcode"] = 0xf029;
	font_awesome_table["question"] = 0xf128;
	font_awesome_table["question-circle"] = 0xf059;
	font_awesome_table["quidditch"] = 0xf458;
	font_awesome_table["quinscape"] = 0xf459;
	font_awesome_table["quora"] = 0xf2c4;
	font_awesome_table["quote-left"] = 0xf10d;
	font_awesome_table["quote-right"] = 0xf10e;
	font_awesome_table["r-project"] = 0xf4f7;
	font_awesome_table["random"] = 0xf074;
	font_awesome_table["ravelry"] = 0xf2d9;
	font_awesome_table["react"] = 0xf41b;
	font_awesome_table["readme"] = 0xf4d5;
	font_awesome_table["rebel"] = 0xf1d0;
	font_awesome_table["receipt"] = 0xf543;
	font_awesome_table["recycle"] = 0xf1b8;
	font_awesome_table["red-river"] = 0xf3e3;
	font_awesome_table["reddit"] = 0xf1a1;
	font_awesome_table["reddit-alien"] = 0xf281;
	font_awesome_table["reddit-square"] = 0xf1a2;
	font_awesome_table["redo"] = 0xf01e;
	font_awesome_table["redo-alt"] = 0xf2f9;
	font_awesome_table["registered"] = 0xf25d;
	font_awesome_table["rendact"] = 0xf3e4;
	font_awesome_table["renren"] = 0xf18b;
	font_awesome_table["reply"] = 0xf3e5;
	font_awesome_table["reply-all"] = 0xf122;
	font_awesome_table["replyd"] = 0xf3e6;
	font_awesome_table["researchgate"] = 0xf4f8;
	font_awesome_table["resolving"] = 0xf3e7;
	font_awesome_table["retweet"] = 0xf079;
	font_awesome_table["rev"] = 0xf5b2;
	font_awesome_table["ribbon"] = 0xf4d6;
	font_awesome_table["road"] = 0xf018;
	font_awesome_table["robot"] = 0xf544;
	font_awesome_table["rocket"] = 0xf135;
	font_awesome_table["rocketchat"] = 0xf3e8;
	font_awesome_table["rockrms"] = 0xf3e9;
	font_awesome_table["route"] = 0xf4d7;
	font_awesome_table["rss"] = 0xf09e;
	font_awesome_table["rss-square"] = 0xf143;
	font_awesome_table["ruble-sign"] = 0xf158;
	font_awesome_table["ruler"] = 0xf545;
	font_awesome_table["ruler-combined"] = 0xf546;
	font_awesome_table["ruler-horizontal"] = 0xf547;
	font_awesome_table["ruler-vertical"] = 0xf548;
	font_awesome_table["rupee-sign"] = 0xf156;
	font_awesome_table["sad-cry"] = 0xf5b3;
	font_awesome_table["sad-tear"] = 0xf5b4;
	font_awesome_table["safari"] = 0xf267;
	font_awesome_table["sass"] = 0xf41e;
	font_awesome_table["save"] = 0xf0c7;
	font_awesome_table["schlix"] = 0xf3ea;
	font_awesome_table["school"] = 0xf549;
	font_awesome_table["screwdriver"] = 0xf54a;
	font_awesome_table["scribd"] = 0xf28a;
	font_awesome_table["search"] = 0xf002;
	font_awesome_table["search-minus"] = 0xf010;
	font_awesome_table["search-plus"] = 0xf00e;
	font_awesome_table["searchengin"] = 0xf3eb;
	font_awesome_table["seedling"] = 0xf4d8;
	font_awesome_table["sellcast"] = 0xf2da;
	font_awesome_table["sellsy"] = 0xf213;
	font_awesome_table["server"] = 0xf233;
	font_awesome_table["servicestack"] = 0xf3ec;
	font_awesome_table["shapes"] = 0xf61f;
	font_awesome_table["share"] = 0xf064;
	font_awesome_table["share-alt"] = 0xf1e0;
	font_awesome_table["share-alt-square"] = 0xf1e1;
	font_awesome_table["share-square"] = 0xf14d;
	font_awesome_table["shekel-sign"] = 0xf20b;
	font_awesome_table["shield-alt"] = 0xf3ed;
	font_awesome_table["ship"] = 0xf21a;
	font_awesome_table["shipping-fast"] = 0xf48b;
	font_awesome_table["shirtsinbulk"] = 0xf214;
	font_awesome_table["shoe-prints"] = 0xf54b;
	font_awesome_table["shopping-bag"] = 0xf290;
	font_awesome_table["shopping-basket"] = 0xf291;
	font_awesome_table["shopping-cart"] = 0xf07a;
	font_awesome_table["shopware"] = 0xf5b5;
	font_awesome_table["shower"] = 0xf2cc;
	font_awesome_table["shuttle-van"] = 0xf5b6;
	font_awesome_table["sign"] = 0xf4d9;
	font_awesome_table["sign-in-alt"] = 0xf2f6;
	font_awesome_table["sign-language"] = 0xf2a7;
	font_awesome_table["sign-out-alt"] = 0xf2f5;
	font_awesome_table["signal"] = 0xf012;
	font_awesome_table["signature"] = 0xf5b7;
	font_awesome_table["simplybuilt"] = 0xf215;
	font_awesome_table["sistrix"] = 0xf3ee;
	font_awesome_table["sitemap"] = 0xf0e8;
	font_awesome_table["sith"] = 0xf512;
	font_awesome_table["skull"] = 0xf54c;
	font_awesome_table["skyatlas"] = 0xf216;
	font_awesome_table["skype"] = 0xf17e;
	font_awesome_table["slack"] = 0xf198;
	font_awesome_table["slack-hash"] = 0xf3ef;
	font_awesome_table["sliders-h"] = 0xf1de;
	font_awesome_table["slideshare"] = 0xf1e7;
	font_awesome_table["smile"] = 0xf118;
	font_awesome_table["smile-beam"] = 0xf5b8;
	font_awesome_table["smile-wink"] = 0xf4da;
	font_awesome_table["smoking"] = 0xf48d;
	font_awesome_table["smoking-ban"] = 0xf54d;
	font_awesome_table["snapchat"] = 0xf2ab;
	font_awesome_table["snapchat-ghost"] = 0xf2ac;
	font_awesome_table["snapchat-square"] = 0xf2ad;
	font_awesome_table["snowflake"] = 0xf2dc;
	font_awesome_table["solar-panel"] = 0xf5ba;
	font_awesome_table["sort"] = 0xf0dc;
	font_awesome_table["sort-alpha-down"] = 0xf15d;
	font_awesome_table["sort-alpha-up"] = 0xf15e;
	font_awesome_table["sort-amount-down"] = 0xf160;
	font_awesome_table["sort-amount-up"] = 0xf161;
	font_awesome_table["sort-down"] = 0xf0dd;
	font_awesome_table["sort-numeric-down"] = 0xf162;
	font_awesome_table["sort-numeric-up"] = 0xf163;
	font_awesome_table["sort-up"] = 0xf0de;
	font_awesome_table["soundcloud"] = 0xf1be;
	font_awesome_table["spa"] = 0xf5bb;
	font_awesome_table["space-shuttle"] = 0xf197;
	font_awesome_table["speakap"] = 0xf3f3;
	font_awesome_table["spinner"] = 0xf110;
	font_awesome_table["splotch"] = 0xf5bc;
	font_awesome_table["spotify"] = 0xf1bc;
	font_awesome_table["spray-can"] = 0xf5bd;
	font_awesome_table["square"] = 0xf0c8;
	font_awesome_table["square-full"] = 0xf45c;
	font_awesome_table["squarespace"] = 0xf5be;
	font_awesome_table["stack-exchange"] = 0xf18d;
	font_awesome_table["stack-overflow"] = 0xf16c;
	font_awesome_table["stamp"] = 0xf5bf;
	font_awesome_table["star"] = 0xf005;
	font_awesome_table["star-half"] = 0xf089;
	font_awesome_table["star-half-alt"] = 0xf5c0;
	font_awesome_table["star-of-life"] = 0xf621;
	font_awesome_table["staylinked"] = 0xf3f5;
	font_awesome_table["steam"] = 0xf1b6;
	font_awesome_table["steam-square"] = 0xf1b7;
	font_awesome_table["steam-symbol"] = 0xf3f6;
	font_awesome_table["step-backward"] = 0xf048;
	font_awesome_table["step-forward"] = 0xf051;
	font_awesome_table["stethoscope"] = 0xf0f1;
	font_awesome_table["sticker-mule"] = 0xf3f7;
	font_awesome_table["sticky-note"] = 0xf249;
	font_awesome_table["stop"] = 0xf04d;
	font_awesome_table["stop-circle"] = 0xf28d;
	font_awesome_table["stopwatch"] = 0xf2f2;
	font_awesome_table["store"] = 0xf54e;
	font_awesome_table["store-alt"] = 0xf54f;
	font_awesome_table["strava"] = 0xf428;
	font_awesome_table["stream"] = 0xf550;
	font_awesome_table["street-view"] = 0xf21d;
	font_awesome_table["strikethrough"] = 0xf0cc;
	font_awesome_table["stripe"] = 0xf429;
	font_awesome_table["stripe-s"] = 0xf42a;
	font_awesome_table["stroopwafel"] = 0xf551;
	font_awesome_table["studiovinari"] = 0xf3f8;
	font_awesome_table["stumbleupon"] = 0xf1a4;
	font_awesome_table["stumbleupon-circle"] = 0xf1a3;
	font_awesome_table["subscript"] = 0xf12c;
	font_awesome_table["subway"] = 0xf239;
	font_awesome_table["suitcase"] = 0xf0f2;
	font_awesome_table["suitcase-rolling"] = 0xf5c1;
	font_awesome_table["sun"] = 0xf185;
	font_awesome_table["superpowers"] = 0xf2dd;
	font_awesome_table["superscript"] = 0xf12b;
	font_awesome_table["supple"] = 0xf3f9;
	font_awesome_table["surprise"] = 0xf5c2;
	font_awesome_table["swatchbook"] = 0xf5c3;
	font_awesome_table["swimmer"] = 0xf5c4;
	font_awesome_table["swimming-pool"] = 0xf5c5;
	font_awesome_table["sync"] = 0xf021;
	font_awesome_table["sync-alt"] = 0xf2f1;
	font_awesome_table["syringe"] = 0xf48e;
	font_awesome_table["font_awesome_table"] = 0xf0ce;
	font_awesome_table["font_awesome_table-tennis"] = 0xf45d;
	font_awesome_table["font_awesome_tablet"] = 0xf10a;
	font_awesome_table["font_awesome_tablet-alt"] = 0xf3fa;
	font_awesome_table["font_awesome_tablets"] = 0xf490;
	font_awesome_table["tachometer-alt"] = 0xf3fd;
	font_awesome_table["tag"] = 0xf02b;
	font_awesome_table["tags"] = 0xf02c;
	font_awesome_table["tape"] = 0xf4db;
	font_awesome_table["tasks"] = 0xf0ae;
	font_awesome_table["taxi"] = 0xf1ba;
	font_awesome_table["teamspeak"] = 0xf4f9;
	font_awesome_table["teeth"] = 0xf62e;
	font_awesome_table["teeth-open"] = 0xf62f;
	font_awesome_table["telegram"] = 0xf2c6;
	font_awesome_table["telegram-plane"] = 0xf3fe;
	font_awesome_table["tencent-weibo"] = 0xf1d5;
	font_awesome_table["terminal"] = 0xf120;
	font_awesome_table["text-height"] = 0xf034;
	font_awesome_table["text-width"] = 0xf035;
	font_awesome_table["th"] = 0xf00a;
	font_awesome_table["th-large"] = 0xf009;
	font_awesome_table["th-list"] = 0xf00b;
	font_awesome_table["theater-masks"] = 0xf630;
	font_awesome_table["themeco"] = 0xf5c6;
	font_awesome_table["themeisle"] = 0xf2b2;
	font_awesome_table["thermometer"] = 0xf491;
	font_awesome_table["thermometer-empty"] = 0xf2cb;
	font_awesome_table["thermometer-full"] = 0xf2c7;
	font_awesome_table["thermometer-half"] = 0xf2c9;
	font_awesome_table["thermometer-quarter"] = 0xf2ca;
	font_awesome_table["thermometer-three-quarters"] = 0xf2c8;
	font_awesome_table["thumbs-down"] = 0xf165;
	font_awesome_table["thumbs-up"] = 0xf164;
	font_awesome_table["thumbtack"] = 0xf08d;
	font_awesome_table["ticket-alt"] = 0xf3ff;
	font_awesome_table["times"] = 0xf00d;
	font_awesome_table["times-circle"] = 0xf057;
	font_awesome_table["tint"] = 0xf043;
	font_awesome_table["tint-slash"] = 0xf5c7;
	font_awesome_table["tired"] = 0xf5c8;
	font_awesome_table["toggle-off"] = 0xf204;
	font_awesome_table["toggle-on"] = 0xf205;
	font_awesome_table["toolbox"] = 0xf552;
	font_awesome_table["tooth"] = 0xf5c9;
	font_awesome_table["trade-federation"] = 0xf513;
	font_awesome_table["trademark"] = 0xf25c;
	font_awesome_table["traffic-light"] = 0xf637;
	font_awesome_table["train"] = 0xf238;
	font_awesome_table["transgender"] = 0xf224;
	font_awesome_table["transgender-alt"] = 0xf225;
	font_awesome_table["trash"] = 0xf1f8;
	font_awesome_table["trash-alt"] = 0xf2ed;
	font_awesome_table["tree"] = 0xf1bb;
	font_awesome_table["trello"] = 0xf181;
	font_awesome_table["tripadvisor"] = 0xf262;
	font_awesome_table["trophy"] = 0xf091;
	font_awesome_table["truck"] = 0xf0d1;
	font_awesome_table["truck-loading"] = 0xf4de;
	font_awesome_table["truck-monster"] = 0xf63b;
	font_awesome_table["truck-moving"] = 0xf4df;
	font_awesome_table["truck-pickup"] = 0xf63c;
	font_awesome_table["tshirt"] = 0xf553;
	font_awesome_table["tty"] = 0xf1e4;
	font_awesome_table["tumblr"] = 0xf173;
	font_awesome_table["tumblr-square"] = 0xf174;
	font_awesome_table["tv"] = 0xf26c;
	font_awesome_table["twitch"] = 0xf1e8;
	font_awesome_table["twitter"] = 0xf099;
	font_awesome_table["twitter-square"] = 0xf081;
	font_awesome_table["typo3"] = 0xf42b;
	font_awesome_table["uber"] = 0xf402;
	font_awesome_table["uikit"] = 0xf403;
	font_awesome_table["umbrella"] = 0xf0e9;
	font_awesome_table["umbrella-beach"] = 0xf5ca;
	font_awesome_table["underline"] = 0xf0cd;
	font_awesome_table["undo"] = 0xf0e2;
	font_awesome_table["undo-alt"] = 0xf2ea;
	font_awesome_table["uniregistry"] = 0xf404;
	font_awesome_table["universal-access"] = 0xf29a;
	font_awesome_table["university"] = 0xf19c;
	font_awesome_table["unlink"] = 0xf127;
	font_awesome_table["unlock"] = 0xf09c;
	font_awesome_table["unlock-alt"] = 0xf13e;
	font_awesome_table["untappd"] = 0xf405;
	font_awesome_table["upload"] = 0xf093;
	font_awesome_table["usb"] = 0xf287;
	font_awesome_table["user"] = 0xf007;
	font_awesome_table["user-alt"] = 0xf406;
	font_awesome_table["user-alt-slash"] = 0xf4fa;
	font_awesome_table["user-astronaut"] = 0xf4fb;
	font_awesome_table["user-check"] = 0xf4fc;
	font_awesome_table["user-circle"] = 0xf2bd;
	font_awesome_table["user-clock"] = 0xf4fd;
	font_awesome_table["user-cog"] = 0xf4fe;
	font_awesome_table["user-edit"] = 0xf4ff;
	font_awesome_table["user-friends"] = 0xf500;
	font_awesome_table["user-graduate"] = 0xf501;
	font_awesome_table["user-lock"] = 0xf502;
	font_awesome_table["user-md"] = 0xf0f0;
	font_awesome_table["user-minus"] = 0xf503;
	font_awesome_table["user-ninja"] = 0xf504;
	font_awesome_table["user-plus"] = 0xf234;
	font_awesome_table["user-secret"] = 0xf21b;
	font_awesome_table["user-shield"] = 0xf505;
	font_awesome_table["user-slash"] = 0xf506;
	font_awesome_table["user-tag"] = 0xf507;
	font_awesome_table["user-tie"] = 0xf508;
	font_awesome_table["user-times"] = 0xf235;
	font_awesome_table["users"] = 0xf0c0;
	font_awesome_table["users-cog"] = 0xf509;
	font_awesome_table["ussunnah"] = 0xf407;
	font_awesome_table["utensil-spoon"] = 0xf2e5;
	font_awesome_table["utensils"] = 0xf2e7;
	font_awesome_table["vaadin"] = 0xf408;
	font_awesome_table["vector-square"] = 0xf5cb;
	font_awesome_table["venus"] = 0xf221;
	font_awesome_table["venus-double"] = 0xf226;
	font_awesome_table["venus-mars"] = 0xf228;
	font_awesome_table["viacoin"] = 0xf237;
	font_awesome_table["viadeo"] = 0xf2a9;
	font_awesome_table["viadeo-square"] = 0xf2aa;
	font_awesome_table["vial"] = 0xf492;
	font_awesome_table["vials"] = 0xf493;
	font_awesome_table["viber"] = 0xf409;
	font_awesome_table["video"] = 0xf03d;
	font_awesome_table["video-slash"] = 0xf4e2;
	font_awesome_table["vimeo"] = 0xf40a;
	font_awesome_table["vimeo-square"] = 0xf194;
	font_awesome_table["vimeo-v"] = 0xf27d;
	font_awesome_table["vine"] = 0xf1ca;
	font_awesome_table["vk"] = 0xf189;
	font_awesome_table["vnv"] = 0xf40b;
	font_awesome_table["volleyball-ball"] = 0xf45f;
	font_awesome_table["volume-down"] = 0xf027;
	font_awesome_table["volume-off"] = 0xf026;
	font_awesome_table["volume-up"] = 0xf028;
	font_awesome_table["vuejs"] = 0xf41f;
	font_awesome_table["walking"] = 0xf554;
	font_awesome_table["wallet"] = 0xf555;
	font_awesome_table["warehouse"] = 0xf494;
	font_awesome_table["weebly"] = 0xf5cc;
	font_awesome_table["weibo"] = 0xf18a;
	font_awesome_table["weight"] = 0xf496;
	font_awesome_table["weight-hanging"] = 0xf5cd;
	font_awesome_table["weixin"] = 0xf1d7;
	font_awesome_table["whatsapp"] = 0xf232;
	font_awesome_table["whatsapp-square"] = 0xf40c;
	font_awesome_table["wheelchair"] = 0xf193;
	font_awesome_table["whmcs"] = 0xf40d;
	font_awesome_table["wifi"] = 0xf1eb;
	font_awesome_table["wikipedia-w"] = 0xf266;
	font_awesome_table["window-close"] = 0xf410;
	font_awesome_table["window-maximize"] = 0xf2d0;
	font_awesome_table["window-minimize"] = 0xf2d1;
	font_awesome_table["window-restore"] = 0xf2d2;
	font_awesome_table["windows"] = 0xf17a;
	font_awesome_table["wine-glass"] = 0xf4e3;
	font_awesome_table["wine-glass-alt"] = 0xf5ce;
	font_awesome_table["wix"] = 0xf5cf;
	font_awesome_table["wolf-pack-battalion"] = 0xf514;
	font_awesome_table["won-sign"] = 0xf159;
	font_awesome_table["wordpress"] = 0xf19a;
	font_awesome_table["wordpress-simple"] = 0xf411;
	font_awesome_table["wpbeginner"] = 0xf297;
	font_awesome_table["wpexplorer"] = 0xf2de;
	font_awesome_table["wpforms"] = 0xf298;
	font_awesome_table["wrench"] = 0xf0ad;
	font_awesome_table["x-ray"] = 0xf497;
	font_awesome_table["xbox"] = 0xf412;
	font_awesome_table["xing"] = 0xf168;
	font_awesome_table["xing-square"] = 0xf169;
	font_awesome_table["y-combinator"] = 0xf23b;
	font_awesome_table["yahoo"] = 0xf19e;
	font_awesome_table["yandex"] = 0xf413;
	font_awesome_table["yandex-international"] = 0xf414;
	font_awesome_table["yelp"] = 0xf1e9;
	font_awesome_table["yen-sign"] = 0xf157;
	font_awesome_table["yoast"] = 0xf2b1;
	font_awesome_table["youtube"] = 0xf167;
	font_awesome_table["youtube-square"] = 0xf431;
	font_awesome_table["zhihu"] = 0xf63f;
    }
    
}
